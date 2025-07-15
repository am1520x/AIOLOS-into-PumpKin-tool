"""
Unit tests for the calculating_rates module.

Tests cover:
- Rate coefficient calculations using modified Arrhenius equation
- Reaction parsing from text files
- Reaction rate calculations from concentrations
"""

import pytest
import numpy as np
import tempfile
import os
import sys
from unittest.mock import mock_open, patch

# Add the Conversion_Tool directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Conversion_Tool'))

from calculating_rates import calc_rate_coeff, parse_reactions, calculate_reaction_rates


class TestCalcRateCoeff:
    """Test cases for the calc_rate_coeff function."""
    
    def test_simple_rate_calculation(self):
        """Test basic rate coefficient calculation."""
        # Simple case: A=1e-10, beta=0, gamma=0, T=300
        result = calc_rate_coeff(1e-10, 0, 0, 300)
        assert abs(result - 1e-10) < 1e-15
        
    def test_temperature_dependence(self):
        """Test temperature dependent rate calculation."""
        A, beta, gamma = 1e-10, 0.5, 1000
        T1, T2 = 300, 600
        
        k1 = calc_rate_coeff(A, beta, gamma, T1)
        k2 = calc_rate_coeff(A, beta, gamma, T2)
        
        # At higher temperature, rate should be higher (for positive gamma)
        assert k2 > k1
        assert k1 > 0
        assert k2 > 0
        
    def test_exponential_term(self):
        """Test the exponential term in the Arrhenius equation."""
        A, beta, gamma, T = 1.0, 0, 1000, 300
        
        expected = 1.0 * np.exp(-1000/300)
        result = calc_rate_coeff(A, beta, gamma, T)
        
        assert abs(result - expected) < 1e-10
        
    def test_power_term(self):
        """Test the temperature power term."""
        A, beta, gamma, T = 1.0, 2.0, 0, 300
        
        expected = 1.0 * (300 ** 2.0)
        result = calc_rate_coeff(A, beta, gamma, T)
        
        assert abs(result - expected) < 1e-10
        
    def test_invalid_inputs(self):
        """Test handling of invalid inputs."""
        # String input should return NaN
        result = calc_rate_coeff("invalid", 0, 0, 300)
        assert np.isnan(result)
        
        # None input should return NaN
        result = calc_rate_coeff(None, 0, 0, 300)
        assert np.isnan(result)
        
    def test_zero_temperature(self):
        """Test behavior with zero temperature."""
        # This should cause division by zero in exp(-gamma/T)
        result = calc_rate_coeff(1e-10, 0, 1000, 0)
        assert np.isnan(result) or np.isinf(result)
        
    def test_negative_temperature(self):
        """Test behavior with negative temperature."""
        result = calc_rate_coeff(1e-10, 0, 1000, -300)
        # Should still compute but may give unexpected physical result
        assert isinstance(result, float)


class TestParseReactions:
    """Test cases for the parse_reactions function."""
    
    def test_parse_single_reaction(self):
        """Test parsing a single reaction."""
        content = "% 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(content)
            temp_file = f.name
            
        try:
            reactions = parse_reactions(temp_file)
            
            assert len(reactions) == 1
            assert reactions[0]['reactants'] == [('H', 1), ('O2', 1)]
            assert reactions[0]['alpha'] == 1.2e-10
            assert reactions[0]['beta'] == 0.5
            assert reactions[0]['gamma'] == 1000
            
        finally:
            os.unlink(temp_file)
            
    def test_parse_multiple_reactions(self):
        """Test parsing multiple reactions."""
        content = """% 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000
% 2 H2 + 1 O -> 1 H2O + 1 H | 3.4e-11 0.0 500
% 1 OH + 1 H -> 1 H2O | 5.6e-12 1.0 200
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(content)
            temp_file = f.name
            
        try:
            reactions = parse_reactions(temp_file)
            
            assert len(reactions) == 3
            
            # Check first reaction
            assert reactions[0]['reactants'] == [('H', 1), ('O2', 1)]
            assert reactions[0]['alpha'] == 1.2e-10
            
            # Check second reaction
            assert reactions[1]['reactants'] == [('H2', 2), ('O', 1)]
            assert reactions[1]['alpha'] == 3.4e-11
            
            # Check third reaction
            assert reactions[2]['reactants'] == [('OH', 1), ('H', 1)]
            assert reactions[2]['alpha'] == 5.6e-12
            
        finally:
            os.unlink(temp_file)
            
    def test_parse_skip_non_reaction_lines(self):
        """Test that non-reaction lines are skipped."""
        content = """# This is a comment
% 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000
This line should be ignored
% 1 OH + 1 H -> 1 H2O | 5.6e-12 1.0 200
Another comment line
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(content)
            temp_file = f.name
            
        try:
            reactions = parse_reactions(temp_file)
            assert len(reactions) == 2
            
        finally:
            os.unlink(temp_file)
            
    def test_file_not_found(self):
        """Test handling of non-existent files."""
        with pytest.raises(FileNotFoundError):
            parse_reactions("non_existent_file.txt")
            
    def test_malformed_reaction_line(self):
        """Test handling of malformed reaction lines."""
        content = "% malformed line without proper format\n"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(content)
            temp_file = f.name
            
        try:
            # This should raise an exception due to malformed format
            with pytest.raises((IndexError, ValueError)):
                parse_reactions(temp_file)
                
        finally:
            os.unlink(temp_file)


class TestCalculateReactionRates:
    """Test cases for the calculate_reaction_rates function."""
    
    def test_single_reaction_rate(self):
        """Test calculation for a single reaction."""
        reactions = [
            {'reactants': [('H', 1), ('O2', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0}
        ]
        densities = {'H': 1e12, 'O2': 1e11}
        
        rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        
        assert len(rates) == 1
        # Expected rate = k * [H] * [O2] = 1e-10 * 1e12 * 1e11 = 1e13
        expected_rate = 1e-10 * 1e12 * 1e11
        assert abs(rates[0] - expected_rate) < 1e-5
        
    def test_multiple_reaction_rates(self):
        """Test calculation for multiple reactions."""
        reactions = [
            {'reactants': [('H', 1), ('O2', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0},
            {'reactants': [('H2', 2)], 'alpha': 2e-11, 'beta': 0, 'gamma': 0},
            {'reactants': [('OH', 1)], 'alpha': 3e-12, 'beta': 0, 'gamma': 0}
        ]
        densities = {'H': 1e12, 'O2': 1e11, 'H2': 1e10, 'OH': 5e9}
        
        rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        
        assert len(rates) == 3
        
        # Check first reaction: rate = 1e-10 * 1e12 * 1e11 = 1e13
        assert abs(rates[0] - 1e13) < 1e5
        
        # Check second reaction: rate = 2e-11 * (1e10)^2 = 2e9
        assert abs(rates[1] - 2e9) < 1e2
        
        # Check third reaction: rate = 3e-12 * 5e9 = 1.5e-2
        assert abs(rates[2] - 1.5e-2) < 1e-5
        
    def test_missing_species_concentration(self):
        """Test handling when species concentration is missing."""
        reactions = [
            {'reactants': [('H', 1), ('UNKNOWN', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0}
        ]
        densities = {'H': 1e12}  # UNKNOWN species not provided
        
        rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        
        assert len(rates) == 1
        # Rate should be zero because UNKNOWN concentration defaults to 0
        assert rates[0] == 0.0
        
    def test_zero_concentrations(self):
        """Test behavior with zero concentrations."""
        reactions = [
            {'reactants': [('H', 1), ('O2', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0}
        ]
        densities = {'H': 0.0, 'O2': 1e11}
        
        rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        
        assert len(rates) == 1
        assert rates[0] == 0.0
        
    def test_temperature_effect_on_rates(self):
        """Test that temperature affects reaction rates."""
        reactions = [
            {'reactants': [('H', 1)], 'alpha': 1e-10, 'beta': 0.5, 'gamma': 1000}
        ]
        densities = {'H': 1e12}
        
        rates_low_T = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        rates_high_T = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=600)
        
        assert len(rates_low_T) == 1
        assert len(rates_high_T) == 1
        # Higher temperature should give higher rate (for positive gamma)
        assert rates_high_T[0] > rates_low_T[0]
        
    def test_stoichiometric_coefficients(self):
        """Test that stoichiometric coefficients are handled correctly."""
        reactions = [
            {'reactants': [('H2', 2), ('O2', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0}
        ]
        densities = {'H2': 1e10, 'O2': 1e11}
        
        rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        
        assert len(rates) == 1
        # Expected rate = k * [H2]^2 * [O2] = 1e-10 * (1e10)^2 * 1e11 = 1e31
        expected_rate = 1e-10 * (1e10)**2 * 1e11
        assert abs(rates[0] - expected_rate) < 1e25


@pytest.mark.integration
class TestIntegratedWorkflow:
    """Integration tests for the complete workflow."""
    
    def test_parse_and_calculate_workflow(self):
        """Test the complete workflow from file parsing to rate calculation."""
        # Create a reaction file
        content = """% 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000
% 2 H2 + 1 O -> 1 H2O + 1 H | 3.4e-11 0.0 500
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(content)
            temp_file = f.name
            
        try:
            # Parse reactions
            reactions = parse_reactions(temp_file)
            
            # Define concentrations
            densities = {'H': 1e12, 'O2': 1e11, 'H2': 1e10, 'O': 1e9}
            
            # Calculate rates
            rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=500)
            
            assert len(reactions) == 2
            assert len(rates) == 2
            assert all(rate >= 0 for rate in rates)
            assert all(isinstance(rate, float) for rate in rates)
            
        finally:
            os.unlink(temp_file)


if __name__ == "__main__":
    pytest.main([__file__])