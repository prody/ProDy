import os
import csv
import pytest

# 1. Keep the skip check at the top
pytest.importorskip("rdkit", allow_module_level=True)

@pytest.fixture
def rdkit_tools():
    """Import RDKit locally so it doesn't crash CI during collection."""
    from rdkit import Chem
    return Chem

@pytest.fixture
def prody_func():
    """Dynamically grab the function from prody."""
    try:
        import prody
        return prody.calcLigandEfficiencyFromSDF
    except (ImportError, AttributeError):
        pytest.fail("prody or calcLigandEfficiencyFromSDF not found.")

@pytest.fixture
def test_sdf(tmp_path, rdkit_tools):
    """Creates a temporary SDF for testing."""
    Chem = rdkit_tools
    sdf_path = tmp_path / "input.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol.SetProp("_Name", "Aspirin_Test")
    mol.SetProp("minimizedAffinity", "-6.5")
    mol.SetProp("minimizedRMSD", "0.5")
    writer.write(mol)
    writer.close()
    
    return str(sdf_path)

## --- Tests ---

def test_prody_import_exists(prody_func):
    """Verify that the tool is accessible from the prody namespace."""
    assert callable(prody_func)

def test_integration_workflow(test_sdf, tmp_path, prody_func):
    """Test the full SDF -> CSV workflow."""
    output_csv = tmp_path / "results.csv"
    
    count = prody_func(
        sdf=test_sdf,
        out=str(output_csv),
        energy_field="minimizedAffinity",
        quiet=True
    )
    
    assert count == 1
    assert os.path.exists(output_csv)
    
    with open(output_csv, 'r') as f:
        reader = csv.DictReader(f)
        row = next(reader)
        assert row['id'] == "Aspirin_Test"
        assert float(row['energy_kcal_per_mol']) == -6.5

def test_limit_logic(test_sdf, tmp_path, prody_func):
    output_csv = tmp_path / "limited.csv"
    count = prody_func(test_sdf, out=str(output_csv), limit=1)
    assert count == 1
