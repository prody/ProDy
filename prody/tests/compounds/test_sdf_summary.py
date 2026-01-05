import os
import csv
import pytest

# 1. Skip all tests if RDKit is not installed
pytest.importorskip("rdkit", allow_module_level=True)

# 2. Import everything from prody
from prody import *
from rdkit import Chem

# 3. Verify the specific function was imported into the namespace
try:
    # This checks if the function exists in the current global namespace
    _ = calcLigandEfficiencyFromSDF
except NameError:
    pytest.fail("calcLigandEfficiencyFromSDF not found in prody namespace. "
                "Ensure the script is correctly integrated into prody.")

@pytest.fixture
def test_sdf(tmp_path):
    """Creates a temporary SDF for testing."""
    sdf_path = tmp_path / "input.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    
    # Create a simple Aspirin molecule
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol.SetProp("_Name", "Aspirin_Test")
    mol.SetProp("minimizedAffinity", "-6.5")
    mol.SetProp("minimizedRMSD", "0.5")
    writer.write(mol)
    writer.close()
    
    return str(sdf_path)

## --- Tests ---

def test_prody_import_exists():
    """Verify that the tool is accessible from the prody namespace."""
    assert callable(calcLigandEfficiencyFromSDF)

def test_integration_workflow(test_sdf, tmp_path):
    """Test the full SDF -> CSV workflow using the prody-imported function."""
    output_csv = tmp_path / "results.csv"
    
    # Run the function imported via 'from prody import *'
    count = calcLigandEfficiencyFromSDF(
        sdf=test_sdf,
        out=str(output_csv),
        energy_field="minimizedAffinity",
        quiet=True
    )
    
    # Assertions
    assert count == 1
    assert os.path.exists(output_csv)
    
    with open(output_csv, 'r') as f:
        reader = csv.DictReader(f)
        row = next(reader)
        assert row['id'] == "Aspirin_Test"
        assert float(row['energy_kcal_per_mol']) == -6.5
        assert float(row['drug_score']) > 0

def test_limit_logic(test_sdf, tmp_path):
    """Verify the limit parameter works when called from prody."""
    # Using the same SDF but asking for a limit
    output_csv = tmp_path / "limited.csv"
    count = calcLigandEfficiencyFromSDF(test_sdf, out=str(output_csv), limit=1)
    assert count == 1
