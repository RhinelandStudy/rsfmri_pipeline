from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    InputMultiPath,
    OutputMultiPath,
    BaseInterface,
    BaseInterfaceInputSpec,
    traits,
    Directory,
    File,
    isdefined,
)
import os
from glob import glob

class FixCleanerInputSpec(CommandLineInputSpec):
    mel_ica = Directory(
        exists=True,
        copyfile=False,
        desc="Melodic output directory or directories",
        argstr="%s",
        position=1,
    )

    trained_wts_file = File(
        exists=True,
        desc="trained-weights (.Rdata) file",
        argstr="%s",
        position=2,
        mandatory=True,
        copyfile=False,
    )

    thresh = traits.Int(
        argstr="%d", desc="Threshold for cleanup.", position=-1, mandatory=True
    )


class FixCleanerOutputSpec(TraitedSpec):

    cleaned_functional_file = File(exists=True, desc="Cleaned session data")
    fix4melview = File(exists=True, desc="fix4melview threshold data")

class FixCleaner(CommandLine):
    """
    Run fix to get clean functional data
    """

    input_spec = FixCleanerInputSpec
    output_spec = FixCleanerOutputSpec
    _cmd = "fix "


    def _list_outputs(self):
        outputs = self.output_spec().get()
        mel_ica = self.inputs.mel_ica
        outputs["cleaned_functional_file"] = os.path.abspath(os.path.join(mel_ica, 'filtered_func_data_clean.nii.gz')) 
        outputs["fix4melview"] = [os.path.abspath(f) for f in glob(os.path.join(mel_ica,'fix4melview*.txt'))][0]
        return outputs

