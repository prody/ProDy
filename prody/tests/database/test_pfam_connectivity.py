import pytest
from prody import searchPfam

def test_searchPfam_connectivity():
    """
    Test connectivity to EBI InterPro API by searching for a 
    standard Pfam accession (PF00001).
    """
    query = "PF00001"
    
    # Perform the search
    results = searchPfam(query, timeout=30)
    
    # 1. Check that we got a result
    assert results is not None, "Pfam search returned None"
    assert isinstance(results, dict), "Results should be a dictionary"
    assert query in results, f"Expected {query} in results"
    
    # 2. Safely check metadata
    metadata = results[query]
    
    # Check for core ID presence
    assert metadata.get("accession") == query, f"Accession mismatch: {metadata.get('accession')}"
    
    # InterPro API 'name' is often inside a 'name' key or a 'description' dict.
    # Let's just verify that 'name' exists somewhere in the metadata keys or values.
    has_name = "name" in metadata or any("name" in str(v).lower() for v in metadata.values())
    assert has_name, "Could not find name information in Pfam metadata"
