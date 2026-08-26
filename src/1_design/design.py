import torch
torch.set_grad_enabled(False)
import argparse

from pytorch_lightning import Trainer

from boltzgen.task.predict.data_from_yaml import FromYamlDataModule, DataConfig
from boltzgen.data.tokenize.tokenizer import Tokenizer
from boltzgen.data.feature.featurizer import Featurizer
from boltzgen.task.predict.writer import DesignWriter
from boltzgen.model.models.boltz import Boltz

MOL_DIR = "/home/anoroozi25/.cache/huggingface/hub/datasets--boltzgen--inference-data/snapshots/c3d36fd276e9caf098c75d4113c6d5eb320b1a4c/mols.zip"
CHECKPOINT = "/home/anoroozi25/.cache/huggingface/hub/models--boltzgen--boltzgen-1/snapshots/c1be29e1f82ffcc72264f64b993c43fb4e0d17f0/boltzgen1_diverse.ckpt"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml_path", type=str, default="boltzgen_inputs/yamls")
    parser.add_argument("--output_dir", type=str, default="boltzgen_outputs/")
    args = parser.parse_args()

    data_config = DataConfig(
        moldir=MOL_DIR,
        multiplicity=100,
        yaml_path=args.yaml_path,
        tokenizer=Tokenizer(
            atomize_modified_residues=False,
            map_to_closest_residue=False
        ),
        featurizer=Featurizer(), # TODO: Add featurizer from config
        backbone_only=False,
        atom14=True,
        atom37=False,
        design=True,
        compute_affinity=False,
        disulfide_prob=1.0,
        disulfide_on=True,
        skip_existing=False,
        skip_offset=0,
        diffusion_samples=1,
        output_dir=args.output_dir
    )
    
    data = FromYamlDataModule(
        data_config,
        batch_size=1,
        num_workers=10,
        pin_memory=True,
        extra_features=None # TODO: Add extra features from config
    )
    
    writer = DesignWriter(
        output_dir=args.output_dir,
        res_atoms_only=False,
        atom14=True,
        atom37=False,
        backbone_only=False,
        write_native=False,
        save_traj=False,
        save_x0_traj=False,
    )
    
    predict_args = {
        "recycling_steps": 3,
        "sampling_steps": 500,
        "diffusion_samples": 10,
        "return_z_feats": False,
        "keys_dict_out": ["s_trunk", "pdistogram", "pbfactor"],
        "confidence_prediction": True,
        "predict_bfactor": False,
    }

    # Load model
    model_module = Boltz.load_from_checkpoint(
        CHECKPOINT,
        strict=True,
        use_ema=False,
        checkpoint_diffusion_conditioning=False,
        map_location="cpu",
        predict_args=predict_args,
        weights_only=False,
    )
    model_module.eval()
    # Print the number of trainable and total model parameters
    total_params = sum(p.numel() for p in model_module.parameters())
    trainable_params = sum(p.numel() for p in model_module.parameters() if p.requires_grad)
    print(f"Total parameters: {total_params}")
    print(f"Trainable parameters: {trainable_params}")

    lightning_trainer = Trainer(
        default_root_dir=args.output_dir,
        strategy="auto",
        callbacks=writer,
        logger=False,
        accelerator="gpu",
        devices=1,
        precision="bf16-mixed",
    )

    # Run training
    lightning_trainer.predict(
        model_module,
        datamodule=data,
    )
