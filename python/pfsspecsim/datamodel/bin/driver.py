import matplotlib.pyplot as plt

from pfs.datamodel.pfsConfig import PfsConfig
from pfs.datamodel.pfsArm import PfsArmSet
from pfs.datamodel.pfsObject import PfsObject, makePfsObject

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def main(pfsConfigId, tract, patch, fiberId=None, dataDir=".", objId=None,
         showPfsArm=False, showPfsArmSet=False, showPfsObject=False):

    pfsConfig = PfsConfig(pfsConfigId, tract, patch)
    pfsConfig.read(dataDir)

    if objId is None:
        if fiberId is None:
            fiberId = 1
    else:
        import numpy as np
        import sys
        try:
            _fiberId = np.where(pfsConfig.objId == objId)[0][0]
        except IndexError:
            print >> sys.stderr, "Unable to find objId %08x in configuration with pfsConfigId 0x%08x" % \
                (objId, pfsConfigId)
            return
        _fiberId += 1                   # 1-indexed
        if fiberId is not None and fiberId != _fiberId:
            print >> sys.stderr, "fiberId %d doesn't match objId %08x's fiber %d" % \
                (fiberId, objId, _fiberId)
            return
        fiberId = _fiberId
        
    objId = pfsConfig.objId[fiberId - 1]

    print "fiberId = %d,  objId = 0x%x" % (fiberId, objId)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    pfsArms = PfsArmSet(visit=1, spectrograph=1, pfsConfigId=pfsConfigId)
    pfsArms.read(dataDir)

    if showPfsArm:
        if True:
            pfsArms.data['r'].plot(fiberId=fiberId)
        else:
            pfsArms.data['r'].plot(fiberId=fiberId,
                                      showFlux=True, showSky=True, showCovar=True, showMask=True)

    if showPfsArmSet:
        if True:
            pfsArms.plot(fiberId=fiberId)
        else:
            pfsArms.plot(fiberId=fiberId,
                         showFlux=True, showSky=True, showCovar=True, showMask=True)

    #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    pfsObject = makePfsObject(objId, [pfsArms])
    pfsObject.write(dataDir)

    npfs = PfsObject(tract, patch, objId, visits=[1])
    npfs.read(dataDir)
    if showPfsObject:
        if True:
            npfs.plot(showFluxTbl=True)
        else:
            npfs.plot(showFlux=True, showSky=True, showCovar=True, showCovar2=True)
            npfs.plot(showFluxTbl=True, showFluxVariance=False)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import argparse
parser = argparse.ArgumentParser(description="Manipulate pfsConfig, pfsArm, and pfsObject files")

parser.add_argument('pfsConfigId', type=str, nargs="?", default="0x788daff233347923",
                    help="Desired pfsConfigId")
parser.add_argument('--dataDir', type=str, default="examples", help="Directory containing data")
parser.add_argument('--fiberId', type=int, default=None, help="Desired fiber")
parser.add_argument('--objId', type=int, default=None, help="Desired object Id")
parser.add_argument('--tract', type=int, default=0, help="Desired tract")
parser.add_argument('--patch', type=str, default="0,0", help="Desired patch")
parser.add_argument('--showPfsArm', action="store_true", help="Plot an pfsArm file?")
parser.add_argument('--showPfsArmSet', action="store_true", help="Plot set of pfsArm files?")
parser.add_argument('--showPfsObject', action="store_true", help="Plot pfsObject file?")

args = parser.parse_args()

main(int(args.pfsConfigId, 16), fiberId=args.fiberId, objId=args.objId,
     tract=args.tract, patch=args.patch, dataDir=args.dataDir,
     showPfsArm=args.showPfsArm, showPfsArmSet=args.showPfsArmSet,
     showPfsObject=args.showPfsObject)
