from sospex.specobj import specCube
import os
import pytest

def checkLoadcube(file):
    path0, file0 = os.path.split(__file__)
    try:
        gcube = specCube(os.path.join(path0,file))
        print(gcube.wcs)
        return True
    except:
        raise NameError('Cube not correctly loaded')     
        return False

def test_loadcube():
    #with pytest.raises(Exception):
    assert checkLoadcube('great.fits')