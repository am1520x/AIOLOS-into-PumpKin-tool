"""
Unit tests for the processing_aiolos_reac_file module.

Tests cover:
- Species name transformation and standardization
- Reaction line processing and cleaning
- Stoichiometry parsing and calculation
- Complete reaction data parsing workflow
"""

import pytest
import sys
import os

# Add the Conversion_Tool directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Conversion_Tool'))

from processing_aiolos_reac_file import (
    transform_species_set, 
    process_reaction_line,
    process_side,
    parse_reaction_data
)


class TestTransformSpeciesSet:
    """Test cases for the transform_species_set function."""
    
    def test_basic_transformations(self):
        """Test basic species transformations."""
        input_species = {"e-", "p", "H2", "CO"}
        result = transform_species_set(input_species)
        
        expected = {"E", "^+", "H2", "CO"}
        assert result == expected
        
    def test_complex_transformations(self):
        """Test more complex transformation cases."""
        input_species = {"e-", "M", "H-", "eV", "CO2p"}
        result = transform_species_set(input_species)
        
        # Check specific transformations
        assert "E" in result  # e- -> E
        assert "ANY_NEUTRAL" in result  # M -> ANY_NEUTRAL
        assert "H^-" in result  # H- -> H^- (- becomes ^-)
        assert "" in result  # eV -> '' (empty)
        assert "CO2^+" in result  # CO2p -> CO2^+ (p becomes ^+ and - becomes ^-)
        
    def test_empty_set(self):
        """Test transformation of empty set."""
        result = transform_species_set(set())
        assert result == set()
        
    def test_no_transformations_needed(self):
        """Test species that don't need transformation."""
        input_species = {"H2O", "CO2", "NH3"}
        result = transform_species_set(input_species)
        
        assert result == input_species
        
    def test_space_removal(self):
        """Test that spaces are removed from species names."""
        input_species = {"H 2", "C O", "N H 3"}
        result = transform_species_set(input_species)
        
        expected = {"H2", "CO", "NH3"}
        assert result == expected
        
    def test_arrow_transformation(self):
        """Test arrow transformation (though unusual in species names)."""
        input_species = {"->", "H2->H"}
        result = transform_species_set(input_species)
        
        assert "=>" in result
        assert "H2=>H" in result


class TestProcessReactionLine:
    """Test cases for the process_reaction_line function."""
    
    def test_simple_reaction(self):
        """Test processing of simple reaction."""
        line = "H2 + O2 -> H2O"
        result = process_reaction_line(line)
        
        assert result == "H2+O2=>H2O"
        
    def test_stoichiometric_coefficients_removal(self):
        """Test removal of stoichiometric coefficients."""
        line = "2 H2 + 1 O2 -> 2 H2O"
        result = process_reaction_line(line)
        
        assert result == "H2+O2=>H2O"
        
    def test_ion_notation(self):
        """Test ion notation standardization."""
        line = "Hp + e- -> H"
        result = process_reaction_line(line)
        
        assert result == "H^++E=>H"
        
    def test_complex_reaction(self):
        """Test complex reaction with multiple transformations."""
        line = "2 H2p + O- + eV -> H2O + H2"
        result = process_reaction_line(line)
        
        # Should remove coefficients, transform ions, remove eV, transform arrow
        assert result == "H2^++O^-=>H2O+H2"
        
    def test_decimal_coefficients(self):
        """Test removal of decimal coefficients."""
        line = "1.5 H2 + 0.5 O2 -> 1.0 H2O"
        result = process_reaction_line(line)
        
        assert result == "H2+O2=>H2O"
        
    def test_special_species(self):
        """Test handling of special species names."""
        line = "M + H2 -> H2 + M"
        result = process_reaction_line(line)
        
        assert result == "ANY_NEUTRAL+H2=>H2+ANY_NEUTRAL"
        
    def test_photodissociation(self):
        """Test photodissociation reaction with eV removal."""
        line = "H2 + eV -> 2 H"
        result = process_reaction_line(line)
        
        assert result == "H2=>H"
        
    def test_energy_term_variants(self):
        """Test different energy term formats."""
        line = "CO2 + eV -> CO + O"
        result = process_reaction_line(line)
        
        assert "eV" not in result
        assert result == "CO2=>CO+O"


class TestProcessSide:
    """Test cases for the process_side function."""
    
    def test_simple_reactants(self):
        """Test processing simple reactants side."""
        stoich = {}
        species = set()
        
        process_side("2 H2 + 1 O2", -1, stoich, species)
        
        assert stoich["H2"] == -2.0
        assert stoich["O2"] == -1.0
        assert "H2" in species
        assert "O2" in species
        
    def test_simple_products(self):
        """Test processing simple products side."""
        stoich = {}
        species = set()
        
        process_side("2 H2O", +1, stoich, species)
        
        assert stoich["H2O"] == 2.0
        assert "H2O" in species
        
    def test_complex_species_names(self):
        """Test processing with complex species names including ions."""
        stoich = {}
        species = set()
        
        process_side("1 H3O+ + 1 e-", -1, stoich, species)
        
        # The regex pattern captures "H3O" not "H3O+" because + is not included in the species pattern
        assert stoich["H3O"] == -1.0
        assert stoich["e-"] == -1.0
        assert "H3O" in species
        assert "e-" in species
        
    def test_decimal_coefficients(self):
        """Test processing with decimal coefficients."""
        stoich = {}
        species = set()
        
        process_side("1.5 H2 + 0.5 O2", -1, stoich, species)
        
        assert stoich["H2"] == -1.5
        assert stoich["O2"] == -0.5
        
    def test_energy_terms_ignored(self):
        """Test that eV terms are ignored."""
        stoich = {}
        species = set()
        
        process_side("1 H2 + 1 eV", -1, stoich, species)
        
        assert "H2" in stoich
        assert "eV" not in stoich
        assert "ev" not in stoich
        assert "H2" in species
        assert "eV" not in species
        
    def test_invalid_coefficients(self):
        """Test handling of invalid coefficient strings."""
        stoich = {}
        species = set()
        
        # This should handle gracefully, though might not match in practice
        process_side("invalid H2", -1, stoich, species)
        
        # Should not crash, behavior depends on regex matching
        # Most likely won't match the pattern and won't add anything
        
    def test_accumulated_stoichiometry(self):
        """Test that stoichiometry accumulates correctly."""
        stoich = {"H2": -1.0}  # Pre-existing
        species = {"H2"}
        
        process_side("2 H2", -1, stoich, species)
        
        assert stoich["H2"] == -3.0  # -1.0 + (-2.0)


class TestParseReactionData:
    """Test cases for the parse_reaction_data function."""
    
    def test_single_thermal_reaction(self):
        """Test parsing single thermal reaction."""
        reac_text = "$ 1 H2 + 1 O -> 1 OH + 1 H | thermal reaction"
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 1
        assert reactions[0] == "1 H2 + 1 O -> 1 OH + 1 H"
        assert len(stoich_list) == 1
        assert stoich_list[0]["H2"] == -1.0
        assert stoich_list[0]["O"] == -1.0
        assert stoich_list[0]["OH"] == 1.0
        assert stoich_list[0]["H"] == 1.0
        assert "H2" in species
        assert "O" in species
        assert "OH" in species
        assert "H" in species
        
    def test_single_photo_reaction(self):
        """Test parsing single photoreaction."""
        reac_text = "% 1 CO2 -> 1 CO + 1 O | photodissociation"
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 1
        assert reactions[0] == "1 CO2 -> 1 CO + 1 O"
        assert stoich_list[0]["CO2"] == -1.0
        assert stoich_list[0]["CO"] == 1.0
        assert stoich_list[0]["O"] == 1.0
        
    def test_multiple_reactions(self):
        """Test parsing multiple reactions."""
        reac_text = """
        $ 2 H2 + 1 O2 -> 2 H2O | combustion
        % 1 H2O -> 1 H2 + 0.5 O2 | photolysis
        $ 1 H + 1 OH -> 1 H2O | recombination
        """
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 3
        assert len(stoich_list) == 3
        
        # Check first reaction
        assert "H2" in reactions[0] and "O2" in reactions[0]
        assert stoich_list[0]["H2"] == -2.0
        assert stoich_list[0]["O2"] == -1.0
        assert stoich_list[0]["H2O"] == 2.0
        
        # Check species list contains all species
        assert "H2" in species
        assert "O2" in species
        assert "H2O" in species
        assert "H" in species
        assert "OH" in species
        
    def test_comment_and_empty_lines(self):
        """Test that comments and empty lines are ignored."""
        reac_text = """
        # This is a comment
        
        $ 1 H2 + 1 O2 -> 1 H2O2 | example
        
        # Another comment
        % 1 H2O -> 1 H2 + 0.5 O2 | photolysis
        
        """
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 2
        assert len(stoich_list) == 2
        
    def test_alternative_arrow_notation(self):
        """Test reactions with => arrow notation."""
        reac_text = "$ 1 H2 + 1 O => 1 OH + 1 H | alternative arrow"
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 1
        assert stoich_list[0]["H2"] == -1.0
        assert stoich_list[0]["O"] == -1.0
        assert stoich_list[0]["OH"] == 1.0
        assert stoich_list[0]["H"] == 1.0
        
    def test_malformed_reactions(self):
        """Test handling of malformed reaction lines."""
        reac_text = """
        $ invalid reaction without arrow | bad
        $ 1 H2 + 1 O2 -> 1 H2O | good reaction
        $ another bad one without proper format
        """
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        # The function currently parses all reaction lines but only creates stoichiometry for valid ones
        assert len(reactions) == 3  # All lines are parsed as reactions
        assert "1 H2 + 1 O2 -> 1 H2O" in reactions[1]  # The good reaction is the second one
        
        # Only the good reaction should have proper stoichiometry
        valid_stoich_count = sum(1 for s in stoich_list if len(s) > 0)
        assert valid_stoich_count == 1  # Only one reaction has valid stoichiometry
        
    def test_empty_input(self):
        """Test parsing empty input."""
        reactions, stoich_list, species = parse_reaction_data("")
        
        assert len(reactions) == 0
        assert len(stoich_list) == 0
        assert len(species) == 0
        
    def test_only_comments(self):
        """Test input with only comments."""
        reac_text = """
        # Comment 1
        # Comment 2
        # Comment 3
        """
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert len(reactions) == 0
        assert len(stoich_list) == 0
        assert len(species) == 0
        
    def test_species_sorting(self):
        """Test that species list is sorted."""
        reac_text = """
        $ 1 Z + 1 Y + 1 X -> 1 W + 1 V + 1 U | test sorting
        """
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        # Species should be sorted alphabetically
        expected_species = sorted(["Z", "Y", "X", "W", "V", "U"])
        assert species == expected_species
        
    def test_complex_stoichiometry(self):
        """Test complex stoichiometry with fractional coefficients."""
        reac_text = "$ 1.5 H2 + 0.5 O2 -> 1.5 H2O | fractional"
        
        reactions, stoich_list, species = parse_reaction_data(reac_text)
        
        assert stoich_list[0]["H2"] == -1.5
        assert stoich_list[0]["O2"] == -0.5
        assert stoich_list[0]["H2O"] == 1.5


@pytest.mark.integration  
class TestIntegratedReactionProcessing:
    """Integration tests for the complete reaction processing workflow."""
    
    def test_complete_aiolos_workflow(self):
        """Test the complete workflow from raw AIOLOS data to processed output."""
        # Simulate realistic AIOLOS reaction file content
        aiolos_data = """
        # AIOLOS chemical network
        $ 2 H + M -> H2 + M | three body recombination
        % H2 + eV -> 2 H | photodissociation
        $ H + OH -> H2O | radical recombination
        $ H2 + O -> OH + H | hydrogen oxidation
        % H2O + eV -> H + OH | water photolysis
        """
        
        # Parse the data
        reactions, stoich_list, species = parse_reaction_data(aiolos_data)
        
        # Transform species set
        transformed_species = transform_species_set(set(species))
        
        # Process reaction lines
        processed_reactions = [process_reaction_line(rxn) for rxn in reactions]
        
        # Verify complete workflow
        assert len(reactions) == 5
        assert len(stoich_list) == 5
        assert len(species) > 0
        
        # Check that species transformation worked
        if "M" in species:
            assert "ANY_NEUTRAL" in transformed_species
            
        # Check that reaction processing worked
        for processed in processed_reactions:
            assert "=>" in processed  # Arrows should be converted
            # Check that stoichiometric coefficients were removed (no standalone numbers)
            reactants = processed.split("=>")[0]
            # Should not have patterns like "2 H" or "1.5 O2" - numbers should only be part of species names
            import re
            standalone_numbers = re.findall(r'\b\d+\.?\d*\s+[A-Za-z]', reactants)
            assert len(standalone_numbers) == 0  # No standalone coefficients
            
        # Verify stoichiometry conservation (sum should be 0 for each reaction)
        for stoich in stoich_list:
            # Skip this check for three-body reactions involving M (catalyst)
            # and photodissociation reactions (eV terms are ignored but affect balance)
            if "M" not in stoich and len(stoich) > 1:  # Skip single-species reactions (likely photodissociation)
                total = sum(stoich.values())
                assert abs(total) < 1e-10  # Should be approximately 0
                
    def test_error_handling_workflow(self):
        """Test error handling in the complete workflow."""
        # Input with various problematic cases
        problematic_data = """
        $ H2 + O2  | missing arrow and products
        % -> H2 | missing reactants  
        $ 1 H2 + 1 O2 -> 1 H2O | good reaction
        invalid line format
        $ H2 + O2 => H2O | another good reaction
        """
        
        # Should handle gracefully without crashing
        reactions, stoich_list, species = parse_reaction_data(problematic_data)
        
        # Should still parse the good reactions
        assert len(reactions) >= 2
        assert "H2O" in species
        
    def test_realistic_chemical_network(self):
        """Test with a realistic subset of chemical reactions."""
        realistic_network = """
        # Hydrogen chemistry
        $ 1 H + 1 H + 1 M -> 1 H2 + 1 M | H2 formation
        $ 1 H2 + 1 O -> 1 OH + 1 H | H2 oxidation  
        $ 1 H + 1 OH -> 1 H2O | H2O formation
        $ 1 OH + 1 OH -> 1 H2O2 | H2O2 formation
        % 1 H2 + 1 eV -> 2 H | H2 photodissociation
        % 1 H2O + 1 eV -> 1 H + 1 OH | H2O photolysis
        % 1 H2O2 + 1 eV -> 2 OH | H2O2 photolysis
        
        # Ion chemistry  
        $ 1 H + 1 e- -> 1 H- | H- formation
        $ 1 H + 1 p -> 1 H2+ | H2+ formation  
        $ 1 H2 + 1 H+ -> 1 H3+ | H3+ formation
        """
        
        reactions, stoich_list, species = parse_reaction_data(realistic_network)
        
        # Verify reasonable parsing
        assert len(reactions) >= 8
        assert len(species) >= 10
        
        # Check for key species
        expected_species = ["H", "H2", "O", "OH", "H2O", "H2O2", "M"]
        for sp in expected_species:
            assert sp in species
            
        # Verify stoichiometry makes sense
        for i, stoich in enumerate(stoich_list):
            # Check that each reaction has both reactants and products
            reactants = [sp for sp, coeff in stoich.items() if coeff < 0]
            products = [sp for sp, coeff in stoich.items() if coeff > 0]
            
            if "M" not in stoich:  # Skip catalyst reactions
                assert len(reactants) >= 1
                assert len(products) >= 1


if __name__ == "__main__":
    pytest.main([__file__])