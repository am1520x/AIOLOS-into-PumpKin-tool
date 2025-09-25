"""
Unit tests for the making_densities_file module.

Tests cover:
- Loading and processing AIOLOS simulation data
- Timestep processing and file output
- Radial profile extraction
- Data validation and error handling
"""

import pytest
import numpy as np
import tempfile
import os
import sys
from unittest.mock import mock_open, patch, MagicMock

# Add the Conversion_Tool directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Conversion_Tool"))

from making_densities_file import (
    get_essential_data,
    process_timesteps,
    process_radial_profile,
)


class TestGetEssentialData:
    """Test cases for the get_essential_data function."""

    def create_mock_data_file(self, n_rows=50):
        """Create mock simulation data array."""
        # Create realistic AIOLOS output data
        # Columns: radius, ?, ?, ?, number_density, ?, ?, ?, ?, ?, ?, temperature, ?
        radius = np.linspace(1.37e8, 5e8, n_rows)  # Planet surface to 5 planet radii
        col1 = np.ones(n_rows) * 0.5
        col2 = np.linspace(1e-3, 1e-6, n_rows)
        col3 = np.ones(n_rows) * 1.0
        number_density = np.logspace(12, 8, n_rows)  # Decreasing with altitude
        col5 = np.ones(n_rows) * 2.0
        col6 = np.ones(n_rows) * 3.0
        col7 = np.ones(n_rows) * 4.0
        col8 = np.ones(n_rows) * 5.0
        col9 = np.ones(n_rows) * 6.0
        col10 = np.ones(n_rows) * 7.0
        temperature = np.linspace(1500, 200, n_rows)  # Decreasing with altitude
        col12 = np.ones(n_rows) * 8.0

        data = np.column_stack(
            [
                radius,
                col1,
                col2,
                col3,
                number_density,
                col5,
                col6,
                col7,
                col8,
                col9,
                col10,
                temperature,
                col12,
            ]
        )
        return data

    @patch("numpy.genfromtxt")
    @patch("os.path.isfile")
    def test_simple_species_list(self, mock_isfile, mock_genfromtxt):
        """Test loading data with simple species name list."""
        # Mock file existence
        mock_isfile.return_value = True

        # Create mock data
        mock_data = self.create_mock_data_file(20)
        mock_genfromtxt.return_value = mock_data

        # Test function
        species = ["H", "He"]
        avg_T, r, all_data = get_essential_data("data/", "sim1", 10, species)

        # Verify results
        assert isinstance(avg_T, np.ndarray)
        assert isinstance(r, np.ndarray)
        assert len(all_data) == 2  # One array per species
        assert len(avg_T) == len(r)
        assert len(avg_T) == 20

        # Check that files were read correctly
        assert mock_genfromtxt.call_count >= 2  # At least first species + others

    @patch("numpy.genfromtxt")
    @patch("os.path.isfile")
    def test_species_data_arrays(self, mock_isfile, mock_genfromtxt):
        """Test loading data with species data arrays (typical format)."""
        mock_isfile.return_value = True
        mock_data = self.create_mock_data_file(15)
        mock_genfromtxt.return_value = mock_data

        # Test with species arrays (format from main pipeline)
        species = [["H", 1, 1, "$\\rm H$"], ["He", 4, 2, "$\\rm He$"]]
        avg_T, r, all_data = get_essential_data(
            "data/", "sim1", 10, species, rplanet=1.37e8
        )

        assert len(all_data) == 2
        assert isinstance(avg_T, np.ndarray)
        assert len(avg_T) == 15

    @patch("numpy.genfromtxt")
    @patch("os.path.isfile")
    def test_planet_radius_normalization(self, mock_isfile, mock_genfromtxt):
        """Test radius normalization with different planet radius settings."""
        mock_isfile.return_value = True
        mock_data = self.create_mock_data_file(10)
        # Set specific planet radius in data
        mock_data[2, 0] = 1.5e8  # Planet radius at index 2
        mock_genfromtxt.return_value = mock_data

        species = ["H"]

        # Test auto-detection (rplanet < 0)
        avg_T, r, all_data = get_essential_data(
            "data/", "sim1", 10, species, rplanet=-1
        )
        assert (
            r[0] == mock_data[0, 0] / 1.5e8
        )  # Should normalize by auto-detected radius

        # Test manual planet radius
        avg_T, r, all_data = get_essential_data(
            "data/", "sim1", 10, species, rplanet=2e8
        )
        assert r[0] == mock_data[0, 0] / 2e8  # Should normalize by specified radius

    @patch("numpy.genfromtxt")
    @patch("os.path.isfile")
    def test_missing_species_file(self, mock_isfile, mock_genfromtxt):
        """Test handling when some species files are missing."""
        # First species exists, second doesn't
        mock_isfile.side_effect = lambda file: "H_t10.dat" in file

        mock_data = self.create_mock_data_file(10)
        mock_genfromtxt.return_value = mock_data

        species = ["H", "MISSING"]
        avg_T, r, all_data = get_essential_data("data/", "sim1", 10, species)

        assert len(all_data) == 2
        # Second species should have zeros (from np.zeros_like)
        assert np.allclose(all_data[1], 0)

    @patch("numpy.genfromtxt")
    def test_file_not_found_error(self, mock_genfromtxt):
        """Test handling when first species file doesn't exist."""
        mock_genfromtxt.side_effect = FileNotFoundError("File not found")

        species = ["H"]
        with pytest.raises(FileNotFoundError):
            get_essential_data("data/", "sim1", 10, species)

    @patch("numpy.genfromtxt")
    @patch("os.path.isfile")
    def test_nan_handling(self, mock_isfile, mock_genfromtxt):
        """Test proper handling of NaN values in data."""
        mock_isfile.return_value = True
        mock_data = self.create_mock_data_file(5)
        # Introduce some NaN values
        mock_data[2, 4] = np.nan  # NaN in number density
        mock_data[3, 11] = np.nan  # NaN in temperature
        mock_genfromtxt.return_value = mock_data

        species = ["H"]
        avg_T, r, all_data = get_essential_data("data/", "sim1", 10, species)

        # Function should handle NaNs gracefully
        assert isinstance(avg_T, np.ndarray)
        assert not np.any(np.isnan(avg_T))  # NaNs should be converted to 0


class TestProcessTimesteps:
    """Test cases for the process_timesteps function."""

    @patch("making_densities_file.get_essential_data")
    def test_basic_timestep_processing(self, mock_get_data):
        """Test basic timestep processing functionality."""
        # Mock return data from get_essential_data
        mock_avg_T = np.array([100, 200, 300, 400, 500])
        mock_r = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
        mock_species_data = [
            np.array(
                [
                    [1e8, 0, 1e12, 0, 100],
                    [1.5e8, 0, 5e11, 0, 200],
                    [2e8, 0, 1e11, 0, 300],
                    [2.5e8, 0, 5e10, 0, 400],
                    [3e8, 0, 1e10, 0, 500],
                ]
            ),  # H data
            np.array(
                [
                    [1e8, 0, 5e11, 0, 100],
                    [1.5e8, 0, 2e11, 0, 200],
                    [2e8, 0, 5e10, 0, 300],
                    [2.5e8, 0, 2e10, 0, 400],
                    [3e8, 0, 5e9, 0, 500],
                ]
            ),  # He data
        ]
        mock_get_data.return_value = (mock_avg_T, mock_r, mock_species_data)

        # Test processing
        timesteps = [0, 5, 10]
        species = [["H", 1, 1, "$\\rm H$"], ["He", 4, 2, "$\\rm He$"]]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "densities.txt")
            output_file_2 = os.path.join(temp_dir, "conditions.txt")

            avg_T_list, data_dict = process_timesteps(
                "data/",
                "sim1",
                timesteps,
                species,
                index=2,
                output_file=output_file,
                output_file_2=output_file_2,
            )

            # Check returned data
            assert len(avg_T_list) == 3  # One per timestep
            assert len(data_dict) == 3  # One per timestep
            assert all(t in data_dict for t in timesteps)
            assert "H" in data_dict[0] and "He" in data_dict[0]

            # Check that files were created
            assert os.path.exists(output_file)
            assert os.path.exists(output_file_2)

    @patch("making_densities_file.get_essential_data")
    def test_timestep_processing_with_errors(self, mock_get_data):
        """Test handling of errors during timestep processing."""

        # First call succeeds, second fails, third succeeds
        def side_effect(*args, **kwargs):
            if kwargs.get("timestep") == 5:
                raise ValueError("Simulated error")
            else:
                mock_avg_T = np.array([100, 200, 300])
                mock_r = np.array([1.0, 1.5, 2.0])
                mock_species_data = [
                    np.array(
                        [
                            [1e8, 0, 1e12, 0, 100],
                            [1.5e8, 0, 5e11, 0, 200],
                            [2e8, 0, 1e11, 0, 300],
                        ]
                    )
                ]
                return mock_avg_T, mock_r, mock_species_data

        mock_get_data.side_effect = side_effect

        timesteps = [0, 5, 10]
        species = [["H", 1, 1, "$\\rm H$"]]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "densities.txt")
            output_file_2 = os.path.join(temp_dir, "conditions.txt")

            avg_T_list, data_dict = process_timesteps(
                "data/",
                "sim1",
                timesteps,
                species,
                output_file=output_file,
                output_file_2=output_file_2,
            )

            # Should have processed all timesteps, with NaN for failed one
            assert len(data_dict) == 3
            assert 5 in data_dict  # Failed timestep should still be in dict
            assert np.isnan(data_dict[5]["H"])  # Should be NaN for failed timestep

    def test_empty_timesteps_list(self):
        """Test behavior with empty timesteps list."""
        species = [["H", 1, 1, "$\\rm H$"]]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "densities.txt")
            output_file_2 = os.path.join(temp_dir, "conditions.txt")

            avg_T_list, data_dict = process_timesteps(
                "data/",
                "sim1",
                [],
                species,
                output_file=output_file,
                output_file_2=output_file_2,
            )

            assert len(avg_T_list) == 0
            assert len(data_dict) == 0


class TestProcessRadialProfile:
    """Test cases for the process_radial_profile function."""

    @patch("making_densities_file.get_essential_data")
    def test_basic_radial_processing(self, mock_get_data):
        """Test basic radial profile processing."""
        # Mock return data
        n_radial = 20
        mock_avg_T = np.linspace(1500, 200, n_radial)  # Temperature profile
        mock_r = np.linspace(1.0, 5.0, n_radial)  # Radius profile
        mock_species_data = [
            np.column_stack(
                [
                    mock_r * 1.37e8,  # radius column
                    np.ones(n_radial),  # dummy column
                    np.logspace(12, 8, n_radial),  # number density (decreasing)
                    np.ones(n_radial),  # dummy column
                    mock_avg_T,  # temperature column
                ]
            ),
            np.column_stack(
                [
                    mock_r * 1.37e8,  # radius column
                    np.ones(n_radial),  # dummy column
                    np.logspace(11, 7, n_radial),  # number density (different species)
                    np.ones(n_radial),  # dummy column
                    mock_avg_T,  # temperature column
                ]
            ),
        ]
        mock_get_data.return_value = (mock_avg_T, mock_r, mock_species_data)

        # Test processing
        species = [["H", 1, 1, "$\\rm H$"], ["He", 4, 2, "$\\rm He$"]]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "radial_densities.txt")
            output_file_2 = os.path.join(temp_dir, "radial_temps.txt")

            avg_T_list, data_dict = process_radial_profile(
                "data/",
                "sim1",
                85,
                species,
                output_file=output_file,
                output_file_2=output_file_2,
            )

            # Check returned data
            assert len(avg_T_list) == n_radial
            assert len(data_dict) == n_radial
            assert all(i in data_dict for i in range(n_radial))

            # Check data structure
            for i in range(n_radial):
                assert "H" in data_dict[i]
                assert "He" in data_dict[i]
                assert isinstance(data_dict[i]["H"], (int, float))
                assert isinstance(data_dict[i]["He"], (int, float))

            # Check that files were created
            assert os.path.exists(output_file)
            assert os.path.exists(output_file_2)

    @patch("making_densities_file.get_essential_data")
    def test_single_species_radial_profile(self, mock_get_data):
        """Test radial profile with single species."""
        n_radial = 5
        mock_avg_T = np.array([1000, 800, 600, 400, 200])
        mock_r = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
        mock_species_data = [
            np.column_stack(
                [
                    mock_r * 1.37e8,
                    np.ones(n_radial),
                    np.array([1e12, 5e11, 1e11, 5e10, 1e10]),  # number densities
                    np.ones(n_radial),
                    mock_avg_T,
                ]
            )
        ]
        mock_get_data.return_value = (mock_avg_T, mock_r, mock_species_data)

        species = ["H"]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "radial_densities.txt")
            output_file_2 = os.path.join(temp_dir, "radial_temps.txt")

            avg_T_list, data_dict = process_radial_profile(
                "data/",
                "sim1",
                85,
                species,
                output_file=output_file,
                output_file_2=output_file_2,
            )

            assert len(avg_T_list) == 5
            assert len(data_dict) == 5
            assert all("H" in data_dict[i] for i in range(5))

    @patch("making_densities_file.get_essential_data")
    def test_file_creation_error(self, mock_get_data):
        """Test handling of file creation errors."""
        mock_avg_T = np.array([100, 200])
        mock_r = np.array([1.0, 2.0])
        mock_species_data = [np.zeros((2, 5))]
        mock_get_data.return_value = (mock_avg_T, mock_r, mock_species_data)

        species = ["H"]

        # Use invalid path to trigger IOError
        invalid_path = "/invalid/path/file.txt"

        with pytest.raises(IOError):
            process_radial_profile(
                "data/",
                "sim1",
                85,
                species,
                output_file=invalid_path,
                output_file_2=invalid_path,
            )


@pytest.mark.integration
class TestIntegratedDensityWorkflow:
    """Integration tests for the complete density processing workflow."""

    @patch("making_densities_file.get_essential_data")
    def test_timesteps_to_radial_consistency(self, mock_get_data):
        """Test that timestep and radial processing give consistent results."""
        # Create consistent mock data
        n_radial = 10
        mock_avg_T = np.linspace(1000, 200, n_radial)
        mock_r = np.linspace(1.0, 3.0, n_radial)
        mock_species_data = [
            np.column_stack(
                [
                    mock_r * 1.37e8,
                    np.ones(n_radial),
                    np.logspace(12, 8, n_radial),
                    np.ones(n_radial),
                    mock_avg_T,
                ]
            )
        ]
        mock_get_data.return_value = (mock_avg_T, mock_r, mock_species_data)

        species = ["H"]
        timestep = 10
        index = 5  # Middle radial bin

        with tempfile.TemporaryDirectory() as temp_dir:
            # Process using timestep method
            ts_file1 = os.path.join(temp_dir, "ts_densities.txt")
            ts_file2 = os.path.join(temp_dir, "ts_temps.txt")
            avg_T_ts, data_ts = process_timesteps(
                "data/",
                "sim1",
                [timestep],
                species,
                index=index,
                output_file=ts_file1,
                output_file_2=ts_file2,
            )

            # Process using radial method
            rad_file1 = os.path.join(temp_dir, "rad_densities.txt")
            rad_file2 = os.path.join(temp_dir, "rad_temps.txt")
            avg_T_rad, data_rad = process_radial_profile(
                "data/",
                "sim1",
                timestep,
                species,
                output_file=rad_file1,
                output_file_2=rad_file2,
            )

            # Results should be consistent
            assert abs(avg_T_ts[0] - avg_T_rad[index]) < 1e-10
            assert abs(data_ts[timestep]["H"] - data_rad[index]["H"]) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__])
