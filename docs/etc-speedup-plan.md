# ETC 高速化の実装計画: ① cos 漸化式 + ② チャンクレベル並列化

## Context

2026-07-10 のプロファイリング調査で、`run_etc` 参照パイプライン
（3.94 s @ n_workers=1 / 1.51 s @ n_workers=3）の実行時間の ~69% が `np.cos`/`np.sin`、
その大半が `psf.spectro_mtf` の Si ディフォーカスループ（(L,1000) 配列への cos 評価 12 回/呼び出し）
であることが判明。BLAS 行列積は ~1% で無関係。scratchpad プロトタイプで以下を実測済み:

- **① cos 加法定理漸化式**（12→6 回の cos）: カーネル 1.46x、偏差 ~9e-16、
  エンドツーエンド全列最大相対偏差 4.6e-14
- **② 行チャンクのスレッド並列**（ufunc は GIL 解放）: 8 分割で 4.9x（効率 61%）

①+② で現行 1.51 s → **0.6–0.9 s** を目標とする。numba（カーネル ~5x）は**保留**
（本計画末尾のメモ参照）。

## 数学的根拠（① は厳密に同一）

ディフォーカス MTF（`psf.py` の Si ディフォーカスループ、gsetc.c:640-661 の忠実な移植）:

$$
\mathrm{MTF}_{\mathrm{defocus}}(\lambda, u) = \frac{\sum_{i=0}^{5} c_i(\lambda) \left[ \frac{1}{6} \cos(f_1 a_i) + \frac{1}{2} \cos(f_2 a_i) + \frac{1}{3} \right]}{\sum_{i=0}^{5} c_i(\lambda)}
$$

$$
a_i(\lambda, u) = \frac{2 \pi (d_0 + i \Delta)}{2 \, n_{\mathrm{Si}}(\lambda) \, F \, p} \, u
, \qquad f_1 = 0.866025404, \quad f_2 = 0.5
$$

（周波数リテラル $f_1, f_2$ は C ソースの値をそのまま使用。 $c_i$ は吸収・厚み因子
`contrib` に、 $\Delta$ は `ddepth` に、 $d_0$ は `_si_defocus_d0` に対応。）

cos の引数は $i$ について**等差**:

$$
f a_i = \alpha_f + i \theta_f, \qquad
\alpha_f = \frac{2 \pi f \, d_0 \, u}{2 n F p}, \qquad
\theta_f = \frac{2 \pi f \, \Delta \, u}{2 n F p}
$$

なので、加法定理の恒等式

$$
\cos(\alpha + (i+1)\theta) = 2 \cos\theta \cos(\alpha + i\theta) - \cos(\alpha + (i-1)\theta)
$$

により、周波数 $f$ ごとに $\cos\alpha_f, \cos(\alpha_f + \theta_f), \cos\theta_f$ の
**3 回の cos 評価**（2 周波数で計 6 回）と安価な乗算・減算だけで $i = 0, \dots, 5$ の全項が得られる。
**実数演算としては恒等変形であり近似はゼロ**。浮動小数点の丸め順序のみが変わり、
実測偏差はカーネル ~9×10⁻¹⁶、e2e ~5×10⁻¹⁴（oracle テストの ~1e-10、
C 参照ゲートの rtol=1.5e-3 に対し余裕 4 桁以上）。
既存の f7728b9 の「exact reassociation (~1e-15)」と同じ性格の変形。

### 検討済み・不採用: Bessel 関数の直接評価

Manual_v5.pdf p.22 のとおり、cos 3 項式は原著者が「速度のために」ジンク関数を近似したもの:

$$
\phi(x) = \frac{J_1(2 \pi x)}{\pi x} \approx \frac{1}{6} \cos(\sqrt{3} \pi x) + \frac{1}{2} \cos(\pi x) + \frac{1}{3}
$$

$\phi$ を `scipy.special.j1` で厳密に評価する案は以下の実測により不採用:

- **速度**: j1 は cos の 2.34 倍/要素（(4096,1000) で 77 ms vs 33 ms）。深度 6 点で
  6 回の j1 評価が必要な上、**Bessel には加法定理に相当する漸化式がない**
  （加法定理は無限級数）ため cos 漸化式のような削減が効かない。
  正味で現行 12-cos よりやや遅く、漸化式版より ~2.3 倍遅い。
- **忠実性**: 近似誤差は $\mathrm{arg} < 2\pi$ で最大 0.012、$\mathrm{arg} > 2\pi$ では
  最大 ~0.99（cos 和は減衰しないがジンクは減衰する）。C 参照出力は cos 近似で
  生成されているため、厳密 Bessel に置き換えると偏差がゲート許容 rtol=1.5e-3 を
  桁違いに超え、**C 参照ゲートと oracle テストが確実に破綻**する。
  CLAUDE.md の「ported quirks を直さない」不変条件に抵触。
- 物理的により正確な MTF が科学的に望ましい場合は、参照出力の再承認を伴う
  上流（ETC 本家）の科学的変更であり、本高速化チケットのスコープ外。

## 実装

### ① `psf.spectro_mtf` の cos 漸化式（`src/pfsspecsim/etc/psf.py`）

Si ディフォーカスの `for i in range(2*n_depth)` ループ（psf.py:161-174）を置換:

- `contrib`/`denom`/`where(depth>thick, 0.3, 1.0)`/`numer/denom` の構造・リテラル
  （`0.1666666667`, `0.5`, `0.3333333333`, `0.866025404`）は**そのまま保持**
- 周波数 `f`（`0.866025404` と `0.5`）ごとに: `theta = f*2π*(ddepth/(nsi*2*fratio*pix))[:,None]*u`、
  `alpha = f*2π*(d0/(nsi*2*fratio*pix))[:,None]*u`、`C0=cos(alpha)`, `C1=cos(alpha+theta)`,
  `T=2*cos(theta)`、以降 `C[i+1]=T*C[i]−C[i−1]` で回して
  `numer += contrib_i[:,None]*(w_f*C_i)` を累積、定数項は `+ 0.3333333333*contrib_i` として従来どおり
- docstring に gsetc.c:640-661 参照と「加法定理による厳密な変形（~1e-15）」の注記を追加
  （既存モジュールの流儀に合わせる）

### ② チャンクレベル並列化（`_parallel.py` + 呼び出し 3 系統）

**新ヘルパー** `src/pfsspecsim/etc/_parallel.py::map_index_chunks(fn, n, chunk_size, n_workers)`:
`range(0, n, chunk_size)` の各チャンクに `fn(start, stop)` を適用。各チャンクは
**互いに素な出力領域に書く**（行スライス／`out[sl]`）ため、`map_arms` と同じ
bit-identity 保証がそのまま成立（順序依存の縮約なし）。`n_workers<=1` または
チャンク数 1 なら従来どおりの逐次ループ（同一コードパス構造）。docstring に
bit-identity の論拠を `map_arms` と同様に明記。

適用箇所（いずれも `params` がスコープ内）:

1. **`engine._map_masked`**（engine.py:153-176、Z_CHUNK=2048 のループ）—
   [OII] カーブ / SNL / MDLF プレスキャン / OII カタログの z スイープを並列化。
   `n_workers` 引数を追加し 4 つの呼び出し元（engine.py:282/360/472/550 相当）から
   `params.n_workers` を渡す。チャンクは `out[sl]` の互いに素な代入 → bit-identical。
2. **`psf.spectro_mtf` に行チャンク並列**: `n_workers: int = 1` kwarg を追加し、
   `L >= 512` かつ `n_workers > 1` のとき行方向を `map_index_chunks` で分割
   （出力 `(L,Nu)` の行は λ ごとに独立 → bit-identical）。
   呼び出し側は当面 **noise.py の 2 箇所のみ**配線
   （noise.py:329 スカイライン、noise.py:378 lam_pix=4096 — noise は
   run_products より前のクリティカルパスなので効果が大きい）。
   snr.py 内部の呼び出しは既存の `mtf=` パススルー構造があるため、
   ①+②-1 の実測後に不足なら追加（YAGNI）。
3. **`EtcParams.n_workers` のデフォルトを 3 → `min(8, os.cpu_count() or 1)`** に変更
   （params.py）。`map_arms` は `min(n_workers, n_arms)` で自然に 3 に飽和するので
   アームレベルの挙動は不変。CLI ヘルプ（cli/etc.py の `--n-workers`）と
   legacy `omp_num_threads` クランプ `[1,3]`（legacy/pfsetc.py）を `[1, 8]` に更新。
   docstring/CLAUDE.md 該当記述の「max 3」も更新。
   オーバーサブスクリプション（アーム 3 × チャンク最大 8）は GIL 解放済み numpy では
   良性と実測済みだが、e2e 計測で悪化が見えたらチャンク側 worker を
   `max(1, n_workers // 2)` に絞る（計測してから決める）。

### テスト（`tests/python/test_psf.py`, `test_engine.py`, `test_parallel` 相当）

- **漸化式 pin テスト**: 従来の 12-cos 直接実装（テスト内に参照実装として転写）と
  `spectro_mtf` を実 Spectrograph（`tests/PFS.20211220.dat`）の Si 両アーム・
  L∈{512, 4096} で比較、rtol=1e-12。既存の「fast-path vs generic-path 等価」
  テストと同じパターン。
- **bit-identity テスト拡張**: 既存の serial/parallel bit-identity テスト
  （test_engine.py）を `n_workers ∈ {1, 2, 8}` グリッドに拡張し、
  チャンク並列込みで全出力テーブルの完全一致を確認。
- **`map_index_chunks` 単体テスト**: 逐次との結果一致、チャンク境界（n がチャンク
  サイズで割り切れる/余る/1 チャンク）のケース。
- 既存の scalar-oracle テスト（~1e-10）はそのまま green になるはず（偏差 ~1e-15）。

## Verification

1. `uv run pytest tests/python -q`（fast suite）
2. `uv run pytest tests/python -m slow -q`（C 参照ゲート — 合格が絶対条件。
   達成済み最大偏差 noise 1.4e-05 / snc 3.2e-04 / snl 4.3e-04 / sno2 7.3e-04 から
   有意に動かないことを確認）
3. ウォールタイム計測: 変更前後で n_workers=1/デフォルトの best-of-3 を比較。
   目標: デフォルト設定で 1.51 s → 0.9 s 以下
4. `etc-port-reviewer` サブエージェントで数値カーネル変更のレビュー
   （CLAUDE.md の運用どおり）

## 保留メモ: numba（③ — 今回は導入しない）

実測: `@njit(cache=True, fastmath=False)` でディフォーカスループ全体
（arg 計算 + 12 cos + 累積）を融合すると **シングルスレッドでもカーネル ~5x**
（L=4096: 220→44 ms）、数値は numpy と**完全一致**（libm 直訳・再結合なし）。
①+② 導入後になお不足する場合の追加カードとして有効。導入時の留意点:

- 初回 JIT コンパイル ~1.7 s（`cache=True` で 2 回目以降は解消）。単発 CLI 実行では
  コンパイルキャッシュが効くまで逆効果になり得る
- プロトタイプで `cache=True` と動的生成クロージャの組で**キャッシュ衝突らしき挙動**を
  観測（parallel/serial 2 変種が同一 qualname を共有した可能性）。導入するなら
  モジュールトップレベルの通常関数として定義すること
- numba↔numpy のバージョン追随ラグ（2026-07 時点では numba 0.66 + numpy 2.4.6 で動作確認済み。
  本リポジトリの `numpy<2.5` ピンとは共存可能）
- 「C コンパイラ不要の純 Python」という v2.0 の設計方針との緊張。導入するなら
  optional dependency（`pip install pfsspecsim[fast]` 等）+ フォールバック実装が筋
- `parallel=True`/`prange` はこの環境ではシングルスレッド njit とほぼ同速だった
  （② のスレッド並列と役割が重複するため、採るなら `parallel=False` で十分）
