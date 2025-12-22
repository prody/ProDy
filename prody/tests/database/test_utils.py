"""Test utilities for database tests with network connectivity checks and fixtures."""

import os
import json
import shutil
from unittest.mock import Mock, patch
from io import BytesIO

from prody import LOGGER

# Global flags to track connectivity
_pfam_connectivity_checked = False
_pfam_is_available = False
_bioexcel_connectivity_checked = False
_bioexcel_is_available = False


def check_pfam_connectivity(timeout=3):
    """
    Check if Pfam/InterPro API is available with a quick connectivity test.
    Returns True if available, False otherwise.
    This should be called once per test session.
    """
    global _pfam_connectivity_checked, _pfam_is_available
    
    if _pfam_connectivity_checked:
        return _pfam_is_available
    
    _pfam_connectivity_checked = True
    
    try:
        import requests
        url = "https://www.ebi.ac.uk/interpro/wwwapi/"
        response = requests.get(url, timeout=timeout)
        _pfam_is_available = response.status_code == 200
        if _pfam_is_available:
            LOGGER.info("Pfam/InterPro API connectivity check: SUCCESS")
        else:
            LOGGER.warn("Pfam/InterPro API connectivity check: FAILED (status {})".format(response.status_code))
    except Exception as e:
        LOGGER.warn("Pfam/InterPro API connectivity check: FAILED ({})".format(str(e)))
        _pfam_is_available = False
    
    return _pfam_is_available


def check_bioexcel_connectivity(timeout=3):
    """
    Check if BioExcel API is available with a quick connectivity test.
    Returns True if available, False otherwise.
    This should be called once per test session.
    """
    global _bioexcel_connectivity_checked, _bioexcel_is_available
    
    if _bioexcel_connectivity_checked:
        return _bioexcel_is_available
    
    _bioexcel_connectivity_checked = True
    
    try:
        import requests
        # Try the mddb-dev API endpoint
        url = "https://irb-dev.mddbr.eu/api/rest/v1/projects/"
        response = requests.head(url, timeout=timeout)
        _bioexcel_is_available = response.status_code in [200, 405]  # 405 means endpoint exists but HEAD not allowed
        if _bioexcel_is_available:
            LOGGER.info("BioExcel API connectivity check: SUCCESS")
        else:
            LOGGER.warn("BioExcel API connectivity check: FAILED (status {})".format(response.status_code))
    except Exception as e:
        LOGGER.warn("BioExcel API connectivity check: FAILED ({})".format(str(e)))
        _bioexcel_is_available = False
    
    return _bioexcel_is_available


def get_fixture_path(fixture_name, subdir=''):
    """Get the full path to a fixture file in the test datafiles directory."""
    import prody
    test_dir = os.path.dirname(prody.tests.__file__)
    datafiles_dir = os.path.join(test_dir, 'datafiles')
    if subdir:
        return os.path.join(datafiles_dir, subdir, fixture_name)
    return os.path.join(datafiles_dir, fixture_name)


def load_pfam_search_fixture(query):
    """Load cached Pfam search results from fixture file."""
    fixture_file = get_fixture_path('{}_search.json'.format(query), 'pfam_fixtures')
    
    if not os.path.exists(fixture_file):
        raise FileNotFoundError("Fixture file not found: {}".format(fixture_file))
    
    with open(fixture_file, 'r') as f:
        data = json.load(f)
    
    return data


def create_mock_parsePDBHeader():
    """
    Create a mock for parsePDBHeader that returns fake polymer data for PDB queries.
    """
    def mock_parse_header(pdb, *keys, **kwargs):
        # Mock polymer data for 6qkc (chain B and I)
        from prody.proteins.header import Polymer, DBRef
        
        if pdb == '6qkc':
            poly_b = Polymer('B')
            dbref_b = DBRef()
            dbref_b.database = 'UniProt'
            dbref_b.accession = 'P19491'  # AMPAR GluA2
            poly_b.dbrefs = [dbref_b]
            
            poly_i = Polymer('I')
            dbref_i = DBRef()
            dbref_i.database = 'UniProt'
            dbref_i.accession = 'Q9JJW0'  # TARP gamma-8
            poly_i.dbrefs = [dbref_i]
            
            if 'polymers' in keys:
                return [poly_b, poly_i]
        
        # Default fallback
        return []
    
    return mock_parse_header


def create_mock_pfam_search(use_fixtures=True, timeout=5):
    """
    Create a mock for searchPfam that uses fixtures.
    
    Args:
        use_fixtures: If True, use cached fixtures. If False, try real network with short timeout.
        timeout: Timeout for network calls if use_fixtures is False.
    
    Returns:
        A function that can replace searchPfam
    """
    def mock_search(query, **kwargs):
        if use_fixtures:
            # Load from fixture
            data = load_pfam_search_fixture(query)
            
            # Process the fixture data the same way searchPfam does
            matches = dict()
            
            if query.startswith('PF'):
                # Pfam ID input
                metadata = data['metadata']
                matches.setdefault(str(query), dict(metadata.items()))
                return matches
            
            # Process results
            for entry in data.get("results", []):
                metadata = entry["metadata"]
                accession = metadata["accession"]
                
                if accession.startswith('PF'):
                    match = matches.setdefault(str(accession), dict(metadata.items()))
                    
                    other_data = entry["proteins"]
                    locations = match.setdefault("locations", [])
                    for item1 in other_data:
                        for key, value in item1.items():
                            if key == "entry_protein_locations":
                                for item2 in value:
                                    for item3 in item2["fragments"]:
                                        new_dict = {}
                                        new_dict["start"] = item3["start"]
                                        new_dict["end"] = item3["end"]
                                        new_dict["score"] = item2.get("score", 1e-10)
                                        locations.append(new_dict)
            
            return matches
        else:
            # Try real network call with timeout
            from prody.database.pfam import searchPfam as real_searchPfam
            kwargs['timeout'] = timeout
            return real_searchPfam(query, **kwargs)
    
    return mock_search


def create_mock_requests_get(use_fixtures=True, timeout=5):
    """
    Create a mock for requests.get that returns fixtures for Pfam/BioExcel.
    This patches at the requests level to catch all network calls.
    """
    import requests
    import gzip
    
    real_get = requests.get
    
    def mock_get(url, **kwargs):
        if use_fixtures:
            # Check if this is a Pfam/InterPro request for search
            if 'interpro/wwwapi/entry' in url and 'annotation=alignment' not in url:
                # Extract query from URL
                if '/protein/uniprot/' in url:
                    query = url.split('/protein/uniprot/')[1].rstrip('/')
                elif '/pfam/' in url:
                    query = url.split('/pfam/')[1].split('?')[0].split('/')[0].rstrip('/')
                else:
                    # Fallback to real request with timeout
                    kwargs['timeout'] = timeout
                    return real_get(url, **kwargs)
                
                try:
                    data = load_pfam_search_fixture(query)
                    
                    # Create a mock response
                    mock_response = Mock()
                    mock_response.content = json.dumps(data).encode('utf-8')
                    mock_response.status_code = 200
                    return mock_response
                except FileNotFoundError:
                    # If fixture not found, try real request with timeout
                    kwargs['timeout'] = timeout
                    return real_get(url, **kwargs)
            
            # Check if this is a Pfam MSA download request
            elif 'interpro/wwwapi/entry' in url and 'annotation=alignment' in url:
                # Extract accession from URL
                parts = url.split('/entry/pfam/')
                if len(parts) > 1:
                    acc = parts[1].split('/')[0]
                    
                    # Determine alignment type from URL
                    if 'annotation=alignment:seed' in url:
                        alignment = 'seed'
                    elif 'annotation=alignment:full' in url:
                        alignment = 'full'
                    else:
                        alignment = 'seed'
                    
                    try:
                        fixture_file = get_fixture_path('{}_{}.sth'.format(acc, alignment), 'pfam_fixtures')
                        if os.path.exists(fixture_file):
                            # Read and gzip the fixture content
                            with open(fixture_file, 'rb') as f:
                                content = f.read()
                            
                            # Gzip it (the real API returns gzipped content)
                            import io
                            buf = io.BytesIO()
                            with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
                                gz.write(content)
                            compressed_content = buf.getvalue()
                            
                            # Create a mock response
                            mock_response = Mock()
                            mock_response.content = compressed_content
                            mock_response.status_code = 200
                            return mock_response
                    except Exception:
                        pass
                
                # Fallback to real request with timeout
                kwargs['timeout'] = timeout
                return real_get(url, **kwargs)
            
            # For non-Pfam requests, use real request with timeout
            kwargs['timeout'] = timeout
            return real_get(url, **kwargs)
        else:
            # Use real requests with timeout
            kwargs['timeout'] = timeout
            return real_get(url, **kwargs)
    
    return mock_get


def create_mock_fetchPfamMSA(use_fixtures=True):
    """Create a mock for fetchPfamMSA that uses fixtures."""
    def mock_fetch(acc, alignment='seed', compressed=False, **kwargs):
        if use_fixtures:
            # Copy fixture to expected location
            fixture_file = get_fixture_path('{}_{}.sth'.format(acc, alignment), 'pfam_fixtures')
            
            if not os.path.exists(fixture_file):
                raise FileNotFoundError("Fixture file not found: {}".format(fixture_file))
            
            folder = kwargs.get('folder', '.')
            outname = kwargs.get('outname', acc)
            from prody.utilities import makePath
            from os.path import join
            
            filepath = join(makePath(folder), outname + '_' + alignment + '.sth')
            
            # Copy the fixture
            shutil.copy(fixture_file, filepath)
            
            from prody.utilities import relpath
            filepath = relpath(filepath)
            LOGGER.info('Pfam MSA for {} is written as {} (from fixture).'.format(acc, filepath))
            
            return filepath
        else:
            # Use real function with short timeout
            from prody.database.pfam import fetchPfamMSA as real_fetchPfamMSA
            kwargs['timeout'] = kwargs.get('timeout', 5)
            return real_fetchPfamMSA(acc, alignment, compressed, **kwargs)
    
    return mock_fetch


def create_mock_ftp_for_pfam_pdbs(use_fixtures=True):
    """
    Create a mock FTP connection for parsePfamPDBs.
    This is more complex as it needs to mock FTP operations.
    """
    if not use_fixtures:
        return None  # Use real FTP
    
    # Mock FTP data - minimal mapping for tests
    mock_pdb_pfam_mapping = """PDB_ID	CHAIN	PDB_START	PDB_END	PFAM_ACC	PFAM_NAME	PFAM_START	PFAM_END
7pj2	A	1	60	PF20446	ATPase_N_2	1	60
7pj2	B	1	60	PF20446	ATPase_N_2	1	60
7pj3	A	1	60	PF20446	ATPase_N_2	1	60
7pj3	B	1	60	PF20446	ATPase_N_2	1	60
7pj4	A	1	60	PF20446	ATPase_N_2	1	60
6yfy	A	264	417	PF01496	V-ATPase_I	1	154
6yfy	A	217	356	PF03223	V-ATPase_H_N	1	140
6yfy	B	264	417	PF01496	V-ATPase_I	1	154
6yfy	B	217	356	PF03223	V-ATPase_H_N	1	140
"""
    
    from ftplib import FTP
    from unittest.mock import MagicMock
    
    class MockFTP:
        def __init__(self, *args, **kwargs):
            pass
        
        def login(self):
            pass
        
        def cwd(self, path):
            pass
        
        def retrbinary(self, cmd, callback):
            # Write the mock data to the callback
            callback(mock_pdb_pfam_mapping.encode('utf-8'))
        
        def quit(self):
            pass
    
    return MockFTP
