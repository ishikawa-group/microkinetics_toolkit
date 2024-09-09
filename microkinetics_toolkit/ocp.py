def set_ocp_calculator():
    """
    Set up a calculator using Neural Network Potential with OCP.
    """
    from fairchem.core.common.relaxation.ase_utils import OCPCalculator
    from fairchem.core.models.model_registry import model_name_to_local_file

    checkpoint_path = model_name_to_local_file("GemNet-dT-S2EFS-OC22", local_cache="./downloaded_checkpoints/")
    calc = OCPCalculator(checkpoint_path=checkpoint_path)

    return calc

