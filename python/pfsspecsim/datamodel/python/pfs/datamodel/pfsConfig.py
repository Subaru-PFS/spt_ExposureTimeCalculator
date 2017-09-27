import os
try:
    import pyfits
except ImportError:
    pyfits = None

from pfs.datamodel.utils import calculate_pfsConfigId, calculate_pfsVisitHash

class PfsConfig(object):
    """A class corresponding to a single pfsConfig file"""
    
    fileNameFormat = "pfsConfig-0x%08x.fits"

    def __init__(self, pfsConfigId=None, tract=None, patch=None,
                 fiberId=None, ra=None, dec=None, catId=None, objId=None,
                 fiberMag=None, mpsCen=None, filterNames=["g", "r", "i", "z", "y"]):
        self.pfsConfigId = pfsConfigId
        self.tract = tract
        self.patch = patch
        
        self.fiberId = fiberId
        self.ra = ra
        self.dec = dec
        self.catId = catId
        self.objId = objId
        self.fiberMag = fiberMag
        self.filterNames = filterNames
        self.mpsCen = mpsCen

        _pfsConfigId = calculate_pfsConfigId(self.fiberId, self.ra, self.dec)

        if self.pfsConfigId is None:
            self.pfsConfigId = _pfsConfigId
        elif _pfsConfigId != 0x0:
            if self.pfsConfigId != _pfsConfigId:
                raise RuntimeError("Mismatch between pfsConfigId == 0x%08x and fiberId/ra/dec -> 0x%08x" %
                                   (self.pfsConfigId, _pfsConfigId))

    def read(self, dirName="."):
        """Read self's pfsConfig file from directory dirName"""
        
        if not pyfits:
            raise RuntimeError("I failed to import pyfits, so cannot read from disk")

        fd = pyfits.open(os.path.join(dirName, self.fileNameFormat % self.pfsConfigId))

        hdu = fd["CONFIG"]
        hdr, data = hdu.header, hdu.data
        self.filterNames = []
        for i in range(5):
            self.filterNames.append(hdr["FILTER%d" % i])
    
        if False:
            for k, v in hdr.items():
                print "%8s %s" % (k, v)
            
        self.fiberId = data['fiberId']
        self.tract = data['tract']
        self.patch = data['patch']
        self.objId = data['objId']
        self.ra = data['ra']
        self.dec = data['dec']
        self.fiberMag = data['fiberMag']
        self.mpsCentroid = data['mps centroid']
        
        assert self.pfsConfigId == calculate_pfsConfigId(self.fiberId, self.ra, self.dec)        

    def write(self, dirName="."):
        if not pyfits:
            raise RuntimeError("I failed to import pyfits, so cannot read from disk")

        for name in ["fiberId", "ra", "dec"]:
            if getattr(self, name) is None:
                raise RuntimeError("I cannot write a pfsConfig file unless %s is provided" % name)

        # even if set in __init__ it might be invalid by now
        _pfsConfigId = calculate_pfsConfigId(self.fiberId, self.ra, self.dec)

        if self.pfsConfigId is None:
            self.pfsConfigId = _pfsConfigId
        else:
            if self.pfsConfigId != _pfsConfigId:
                raise RuntimeError("Mismatch between pfsConfigId == 0x%08x and fiberId/ra/dec -> 0x%08x" %
                                   (self.pfsConfigId, _pfsConfigId))

        hdus = pyfits.HDUList()

        hdr = pyfits.Header()
        hdu = pyfits.PrimaryHDU(header=hdr)
        hdr.update()
        hdus.append(hdu)

        # catId, objId, ra, dec, fiber flux, MPS centroid
        hdr = pyfits.Header()
        for i, b in enumerate(self.filterNames):
            hdr["FILTER%d" % i] = b
        hdr.update(INHERIT=True)
        
        hdu = pyfits.BinTableHDU.from_columns([
            pyfits.Column(name = 'fiberId', format = 'J', array=self.fiberId),
            pyfits.Column(name = 'catId', format = 'J', array=self.catId),
            pyfits.Column(name = 'tract', format = 'J', array=self.tract),
            pyfits.Column(name = 'patch', format = 'A3', array=self.patch),
            pyfits.Column(name = 'objId', format = 'K', array=self.objId),
            pyfits.Column(name = 'ra', format = 'E', array=self.ra),
            pyfits.Column(name = 'dec', format = 'E', array=self.dec),
            pyfits.Column(name = 'fiberMag', format = '5E', array=self.fiberMag),
            pyfits.Column(name = 'MPS centroid', format = '2E', array=self.mpsCen)
        ], hdr)
        hdu.name = 'CONFIG'
        hdus.append(hdu)

        # clobber=True in writeto prints a message, so use open instead
        fileName = self.fileNameFormat % (self.pfsConfigId)
        with open(os.path.join(dirName, fileName), "w") as fd:
            hdus.writeto(fd)            
