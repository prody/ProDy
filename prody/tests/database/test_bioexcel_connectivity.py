import pytest
import os
from prody import fetchBioexcelPDB

def is_server_up():
    try:
        # Quick 3-second check to see if the API base is responding
        response = requests.head("https://bioexcel-cv19.bsc.es/api/rest/v1/", timeout=3)
        return response.status_code == 200
    except:
        return False

@pytest.mark.skipif(not is_server_up(), reason="BioExcel server is unreachable")
def test_fetchBioexcelPDB_connectivity(tmpdir):
    """
    Test connectivity to BioExcel-CV19 / MDDB API.
    """
    acc = "MCV1900001"
    temp_folder = str(tmpdir)
    
    # We do NOT pass timeout here because the source code 
    # has a bug that crashes if timeout is passed as a keyword.
    # It will use the internal default of 200 seconds.
    try:
       # filepath = fetchBioexcelPDB(acc, folder=temp_folder)
        filepath = fetchBioexcelPDB(acc, folder=temp_folder, db='cv19', selection='backbone')        
        assert os.path.exists(filepath), f"File not found at {filepath}"
        assert os.path.getsize(filepath) > 0, "Downloaded file is empty"
        
        with open(filepath, 'r') as f:
            content = f.read(100)
            assert any(tag in content for tag in ["HEADER", "TITLE", "REMARK", "ATOM"]), \
                "File content does not look like a PDB."

    except Exception as e:
        raise e
