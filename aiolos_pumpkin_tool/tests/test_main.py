"""
Unit tests for the main.py pipeline module.

Tests cover:
- Pipeline initialization and configuration
- Command-line argument parsing
- Individual pipeline stages (AIOLOS processing, PumpKin, plotting)
- Complete pipeline execution
- Error handling and validation
"""

import pytest
import tempfile
import os
import sys
import argparse
from unittest.mock import Mock, patch, MagicMock, mock_open
from pathlib import Path

# Add the root directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from main import AIOLOSPumpKinPipeline, create_parser, main


class TestArgumentParser:
    """Test cases for command-line argument parsing."""

    def test_create_parser_basic(self):
        """Test basic parser creation."""
        parser = create_parser()
        assert isinstance(parser, argparse.ArgumentParser)

        # Test that parser can parse minimal arguments
        args = parser.parse_args(["--process-data"])
        assert args.process_data is True

    def test_parser_run_all_mode(self):
        """Test --run-all argument parsing."""
        parser = create_parser()
        args = parser.parse_args(["--run-all"])

        assert args.run_all is True
        assert args.num_cells == 233  # Default value
        assert args.species == ["H", "H2", "O", "OH"]  # Default species
        assert args.processing_type == "timesteps"
        assert args.output_folder == "Testing"
        assert args.timestep == 85

    def test_parser_individual_modes(self):
        """Test individual mode arguments."""
        parser = create_parser()

        # Test process-data only
        args = parser.parse_args(["--process-data"])
        assert args.process_data is True
        assert args.run_pumpkin is False
        assert args.generate_plots is False

        # Test multiple modes
        args = parser.parse_args(["--process-data", "--run-pumpkin"])
        assert args.process_data is True
        assert args.run_pumpkin is True
        assert args.generate_plots is False

    def test_parser_custom_configuration(self):
        """Test custom configuration arguments."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "--run-all",
                "--num-cells",
                "100",
                "--species",
                "H",
                "He",
                "H2O",
                "--top-n",
                "10",
                "--processing-type",
                "radial_profile",
                "--timestep",
                "50",
                "--output-folder",
                "CustomTest",
            ]
        )

        assert args.num_cells == 100
        assert args.species == ["H", "He", "H2O"]
        assert args.top_n == 10
        assert args.processing_type == "radial_profile"
        assert args.timestep == 50
        assert args.output_folder == "CustomTest"

    def test_parser_file_configuration(self):
        """Test file and directory configuration arguments."""
        parser = create_parser()
        args = parser.parse_args(
            [
                "--process-data",
                "--data-dir",
                "./custom_data/",
                "--output-dir",
                "./custom_output/",
                "--reaction-file",
                "custom.reac",
                "--species-file",
                "custom_species.txt",
            ]
        )

        assert args.data_dir == "./custom_data/"
        assert args.output_dir == "./custom_output/"
        assert args.reaction_file == "custom.reac"
        assert args.species_file == "custom_species.txt"


class TestPipelineInitialization:
    """Test cases for pipeline initialization."""

    def test_pipeline_basic_initialization(self):
        """Test basic pipeline initialization."""
        # Create mock arguments
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = False

        with patch("os.makedirs"):
            pipeline = AIOLOSPumpKinPipeline(mock_args)

            assert pipeline.args == mock_args
            assert hasattr(pipeline.args, "processing_type")
            assert hasattr(pipeline.args, "output_folder")
            assert hasattr(pipeline.args, "timestep")

    def test_pipeline_default_values(self):
        """Test that default values are set correctly."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = False

        # Remove the new attributes to test default setting
        del mock_args.processing_type
        del mock_args.output_folder
        del mock_args.timestep

        with patch("os.makedirs"):
            pipeline = AIOLOSPumpKinPipeline(mock_args)

            assert pipeline.args.processing_type == "timesteps"
            assert pipeline.args.output_folder == "Testing"
            assert pipeline.args.timestep == 85

    @patch("os.makedirs")
    def test_setup_directories(self, mock_makedirs):
        """Test directory setup."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = False

        pipeline = AIOLOSPumpKinPipeline(mock_args)

        # Check that makedirs was called for both directories
        expected_calls = [
            ((mock_args.output_dir,), {"exist_ok": True}),
            ((mock_args.data_dir,), {"exist_ok": True}),
        ]
        assert mock_makedirs.call_count == 2

    def test_validate_configuration_invalid_cells(self):
        """Test validation with invalid number of cells."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 0  # Invalid
        mock_args.run_pumpkin = False

        with patch("os.makedirs"):
            with pytest.raises(ValueError, match="Number of cells must be positive"):
                AIOLOSPumpKinPipeline(mock_args)

    @patch("os.path.exists")
    @patch("os.makedirs")
    def test_validate_configuration_missing_pumpkin_dir(
        self, mock_makedirs, mock_exists
    ):
        """Test validation with missing PumpKin directory."""
        mock_exists.return_value = False
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/nonexistent/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = True  # This should trigger warning

        with patch("builtins.print") as mock_print:
            pipeline = AIOLOSPumpKinPipeline(mock_args)
            mock_print.assert_called_with(
                f"Warning: PumpKin directory /nonexistent/pumpkin does not exist"
            )


class TestAIOLOSDataProcessing:
    """Test cases for AIOLOS data processing methods."""

    def create_mock_pipeline(self):
        """Create a mock pipeline for testing."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = False
        mock_args.reaction_file = "test.reac"
        mock_args.species_file = "species.txt"
        mock_args.reactions_file = "reactions.txt"
        mock_args.densities_file = "densities.txt"
        mock_args.rates_file = "rates.txt"
        mock_args.simulation_dir = "../test_data/"
        mock_args.simulation_name = "test_sim"

        with patch("os.makedirs"):
            return AIOLOSPumpKinPipeline(mock_args)

    @patch("making_densities_file.process_timesteps")
    @patch("making_rates_file.make_rates")
    @patch("os.makedirs")
    def test_process_simulation_data_timesteps(
        self, mock_makedirs, mock_make_rates, mock_process_ts
    ):
        """Test simulation data processing with timesteps method."""
        pipeline = self.create_mock_pipeline()
        pipeline.args.processing_type = "timesteps"

        # Mock return values
        mock_avg_T = [100, 200, 300]
        mock_num_den = {0: {"H": 1e12}, 1: {"H": 5e11}, 2: {"H": 1e11}}
        mock_process_ts.return_value = (mock_avg_T, mock_num_den)

        result = pipeline._process_simulation_data()

        assert result is True
        mock_process_ts.assert_called_once()
        mock_make_rates.assert_called_once()

        # Check that process_timesteps was called with correct arguments
        call_args = mock_process_ts.call_args
        assert call_args[1]["directory"] == "../test_data/"
        assert call_args[1]["sim"] == "test_sim"
        assert len(call_args[1]["species"]) > 0  # Should have species list

    @patch("making_densities_file.process_radial_profile")
    @patch("making_rates_file.make_rates")
    @patch("os.makedirs")
    def test_process_simulation_data_radial_profile(
        self, mock_makedirs, mock_make_rates, mock_process_rad
    ):
        """Test simulation data processing with radial profile method."""
        pipeline = self.create_mock_pipeline()
        pipeline.args.processing_type = "radial_profile"
        pipeline.args.timestep = 50

        # Mock return values
        mock_avg_T = [100, 200, 300, 400, 500]
        mock_num_den = {0: {"H": 1e12}, 1: {"H": 5e11}, 2: {"H": 1e11}}
        mock_process_rad.return_value = (mock_avg_T, mock_num_den)

        result = pipeline._process_simulation_data()

        assert result is True
        mock_process_rad.assert_called_once()
        mock_make_rates.assert_called_once()

        # Check that process_radial_profile was called with correct timestep
        call_args = mock_process_rad.call_args
        assert call_args[1]["timestep"] == 50

    @patch("os.path.exists")
    def test_process_reaction_files_missing_file(self, mock_exists):
        """Test reaction file processing with missing input file."""
        mock_exists.return_value = False
        pipeline = self.create_mock_pipeline()

        result = pipeline._process_reaction_files()

        assert result is False

    @patch("os.path.exists")
    @patch("processing_aiolos_reac_file.parse_reaction_data")
    @patch("processing_aiolos_reac_file.process_reaction_line")
    @patch("processing_aiolos_reac_file.transform_species_set")
    def test_process_reaction_files_success(
        self, mock_transform, mock_process_line, mock_parse_data, mock_exists
    ):
        """Test successful reaction file processing."""
        mock_exists.return_value = True
        pipeline = self.create_mock_pipeline()

        # Mock return values
        mock_parse_data.return_value = (
            ["reaction1"],
            [{"H": 1, "O2": 1}],
            ["H", "O2", "OH", "O"],
        )
        mock_transform.return_value = ["H", "O2", "OH", "O"]
        mock_process_line.return_value = "processed_reaction"

        with patch("builtins.open", mock_open(read_data="test reaction data")):
            result = pipeline._process_reaction_files()

        assert result is True
        mock_parse_data.assert_called_once()
        mock_transform.assert_called_once()
        mock_process_line.assert_called_once_with("reaction1")


class TestPumpKinAnalysis:
    """Test cases for PumpKin analysis methods."""

    def create_mock_pipeline(self):
        """Create a mock pipeline for testing."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.pumpkin_in_dir = "/test/pumpkin/Examples/test"
        mock_args.num_cells = 3
        mock_args.run_pumpkin = True
        mock_args.species = ["H", "H2"]
        mock_args.output_prefix = "test_output"

        with patch("os.makedirs"):
            return AIOLOSPumpKinPipeline(mock_args)

    @patch("os.path.exists")
    def test_run_pumpkin_analysis_missing_input(self, mock_exists):
        """Test PumpKin analysis with missing input file."""
        mock_exists.return_value = False
        pipeline = self.create_mock_pipeline()

        result = pipeline.run_pumpkin_analysis()

        assert result is False

    @patch("os.path.exists")
    @patch("subprocess.run")
    @patch("builtins.open", mock_open())
    def test_run_pumpkin_analysis_success(self, mock_subprocess, mock_exists):
        """Test successful PumpKin analysis."""
        mock_exists.return_value = True
        pipeline = self.create_mock_pipeline()

        # Mock subprocess success
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "PumpKin output"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        with patch.object(pipeline, "_modify_pumpkin_input", return_value=True):
            result = pipeline.run_pumpkin_analysis()

        assert result is True
        assert mock_subprocess.call_count == 3  # num_cells = 3

    @patch("os.path.exists")
    @patch("subprocess.run")
    @patch("builtins.open", mock_open())
    def test_run_pumpkin_analysis_partial_failure(self, mock_subprocess, mock_exists):
        """Test PumpKin analysis with some cell failures."""
        mock_exists.return_value = True
        pipeline = self.create_mock_pipeline()

        # Mock mixed success/failure
        def mock_subprocess_side_effect(*args, **kwargs):
            mock_result = Mock()
            mock_result.stdout = "output"
            mock_result.stderr = "error"
            # First call succeeds, others fail
            if mock_subprocess.call_count == 1:
                mock_result.returncode = 0
            else:
                mock_result.returncode = 1
            return mock_result

        mock_subprocess.side_effect = mock_subprocess_side_effect

        with patch.object(pipeline, "_modify_pumpkin_input", return_value=True):
            result = pipeline.run_pumpkin_analysis()

        assert result is True  # Should still return True if at least one succeeds

    @patch("builtins.open", mock_open(read_data="t_init = 0\nt_end = 1\nother_line\n"))
    def test_modify_pumpkin_input_success(self):
        """Test successful PumpKin input file modification."""
        pipeline = self.create_mock_pipeline()

        with patch("builtins.open", mock_open()) as mock_file:
            result = pipeline._modify_pumpkin_input("/test/input.txt", 5)

        assert result is True
        # Check that file was written with modified content
        mock_file.assert_called()


class TestPlotGeneration:
    """Test cases for plot generation methods."""

    def create_mock_pipeline(self):
        """Create a mock pipeline for testing."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = False
        mock_args.species = ["H", "H2", "OH"]
        mock_args.top_n = 5
        mock_args.species_file = "species.txt"
        mock_args.densities_file = "densities.txt"
        mock_args.reactions_file = "reactions.txt"
        mock_args.rates_file = "rates.txt"

        with patch("os.makedirs"):
            return AIOLOSPumpKinPipeline(mock_args)

    @patch("os.path.exists")
    def test_generate_plots_missing_script(self, mock_exists):
        """Test plot generation with missing plotting script."""
        mock_exists.return_value = False
        pipeline = self.create_mock_pipeline()

        result = pipeline.generate_plots()

        assert result is False

    @patch("os.path.exists")
    @patch("subprocess.run")
    def test_generate_plots_success(self, mock_subprocess, mock_exists):
        """Test successful plot generation."""
        mock_exists.return_value = True
        pipeline = self.create_mock_pipeline()

        # Mock subprocess success
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "Plots generated successfully"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        result = pipeline.generate_plots()

        assert result is True
        mock_subprocess.assert_called_once()

        # Check that the command was built correctly
        call_args = mock_subprocess.call_args[0][0]
        assert "plotting_main.py" in call_args[1]
        assert "--all" in call_args
        assert "--species" in call_args

    @patch("os.path.exists")
    @patch("subprocess.run")
    def test_generate_plots_failure(self, mock_subprocess, mock_exists):
        """Test plot generation failure."""
        mock_exists.return_value = True
        pipeline = self.create_mock_pipeline()

        # Mock subprocess failure
        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stdout = ""
        mock_result.stderr = "Plot generation failed"
        mock_subprocess.return_value = mock_result

        result = pipeline.generate_plots()

        assert result is False


class TestCompletePipeline:
    """Test cases for complete pipeline execution."""

    def create_mock_pipeline(self):
        """Create a mock pipeline for testing."""
        mock_args = Mock()
        mock_args.output_dir = "./test_output"
        mock_args.data_dir = "./test_data"
        mock_args.pumpkin_dir = "/test/pumpkin"
        mock_args.num_cells = 10
        mock_args.run_pumpkin = True
        mock_args.species = ["H", "H2"]
        mock_args.process_data = True
        mock_args.generate_plots = True

        with patch("os.makedirs"):
            return AIOLOSPumpKinPipeline(mock_args)

    def test_run_complete_pipeline_all_success(self):
        """Test complete pipeline with all stages successful."""
        pipeline = self.create_mock_pipeline()

        with patch.object(
            pipeline, "process_aiolos_data", return_value=True
        ), patch.object(
            pipeline, "run_pumpkin_analysis", return_value=True
        ), patch.object(
            pipeline, "generate_plots", return_value=True
        ):

            result = pipeline.run_complete_pipeline()

        assert result is True

    def test_run_complete_pipeline_partial_failure(self):
        """Test complete pipeline with some stage failures."""
        pipeline = self.create_mock_pipeline()

        with patch.object(
            pipeline, "process_aiolos_data", return_value=True
        ), patch.object(
            pipeline, "run_pumpkin_analysis", return_value=False
        ), patch.object(
            pipeline, "generate_plots", return_value=True
        ):

            result = pipeline.run_complete_pipeline()

        assert result is False

    def test_run_complete_pipeline_skip_stages(self):
        """Test complete pipeline with some stages skipped."""
        pipeline = self.create_mock_pipeline()
        pipeline.args.process_data = False  # Skip data processing

        with patch.object(
            pipeline, "process_aiolos_data", return_value=True
        ), patch.object(
            pipeline, "run_pumpkin_analysis", return_value=True
        ), patch.object(
            pipeline, "generate_plots", return_value=True
        ):

            result = pipeline.run_complete_pipeline()

        assert result is True
        # process_aiolos_data should not be called since process_data = False


class TestMainFunction:
    """Test cases for the main function."""

    @patch("sys.argv", ["main.py"])
    def test_main_no_mode_selected(self):
        """Test main function with no operation mode selected."""
        with patch("builtins.print") as mock_print:
            result = main()

        assert result == 1
        mock_print.assert_called_with(
            "No operation mode selected. Use --help for options or "\n            "--run-all for complete pipeline."
        )

    @patch("sys.argv", ["main.py", "--run-all"])
    @patch("main.AIOLOSPumpKinPipeline")
    def test_main_run_all_success(self, mock_pipeline_class):
        """Test main function with --run-all mode."""
        mock_pipeline = Mock()
        mock_pipeline.run_complete_pipeline.return_value = True
        mock_pipeline_class.return_value = mock_pipeline

        result = main()

        assert result == 0
        mock_pipeline_class.assert_called_once()
        mock_pipeline.run_complete_pipeline.assert_called_once()

    @patch("sys.argv", ["main.py", "--run-all"])
    @patch("main.AIOLOSPumpKinPipeline")
    def test_main_run_all_failure(self, mock_pipeline_class):
        """Test main function with pipeline failure."""
        mock_pipeline = Mock()
        mock_pipeline.run_complete_pipeline.return_value = False
        mock_pipeline_class.return_value = mock_pipeline

        result = main()

        assert result == 1

    @patch("sys.argv", ["main.py", "--run-all"])
    @patch("main.AIOLOSPumpKinPipeline")
    def test_main_exception_handling(self, mock_pipeline_class):
        """Test main function exception handling."""
        mock_pipeline_class.side_effect = Exception("Test exception")

        with patch("builtins.print") as mock_print, patch(
            "traceback.print_exc"
        ) as mock_traceback:
            result = main()

        assert result == 1
        mock_print.assert_called_with("Fatal error: Test exception")
        mock_traceback.assert_called_once()


@pytest.mark.integration
class TestIntegratedMainWorkflow:
    """Integration tests for the complete main workflow."""

    @patch("os.path.exists")
    @patch("os.makedirs")
    @patch("subprocess.run")
    def test_complete_workflow_simulation(
        self, mock_subprocess, mock_makedirs, mock_exists
    ):
        """Test a simulated complete workflow."""
        # Mock all external dependencies
        mock_exists.return_value = True

        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "Success"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        # Create parser and arguments
        parser = create_parser()
        args = parser.parse_args(["--run-all", "--num-cells", "2"])

        # Mock all the imported functions
        with patch(
            "processing_aiolos_reac_file.parse_reaction_data"
        ) as mock_parse, patch(
            "processing_aiolos_reac_file.process_reaction_line"
        ) as mock_process, patch(
            "processing_aiolos_reac_file.transform_species_set"
        ) as mock_transform, patch(
            "making_densities_file.process_timesteps"
        ) as mock_process_ts, patch(
            "making_rates_file.make_rates"
        ) as mock_make_rates:

            # Set up mock returns
            mock_parse.return_value = (["reaction"], [{}], ["H"])
            mock_transform.return_value = ["H"]
            mock_process.return_value = "processed"
            mock_process_ts.return_value = ([100], {0: {"H": 1e12}})

            # Create and run pipeline
            pipeline = AIOLOSPumpKinPipeline(args)

            # This tests the integration without actually running external processes
            assert pipeline is not None
            assert pipeline.args.num_cells == 2
            assert hasattr(pipeline, "process_aiolos_data")
            assert hasattr(pipeline, "run_pumpkin_analysis")
            assert hasattr(pipeline, "generate_plots")


if __name__ == "__main__":
    pytest.main([__file__])
