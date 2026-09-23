# Main Simulation

Role in paper:

- Original submitted simulation comparing Raw, ComBat, CovBat, and GG-ComBat.
- Evaluates cross-site mean, variance, and covariance alignment.
- Evaluates residual site prediction and leave-one-site-out biological/sex prediction.

Authoritative source snapshot:

- `original/harm_ml.Rmd`

Submitted result files:

- `results/submitted/distances.csv`
- `results/submitted/site_pred.csv`
- `results/submitted/sex_pred_loso.csv`
- `results/submitted/gamma_long.csv`
- `results/submitted/delta2_long.csv`
- `results/submitted/run_meta.csv`

Important revision note:

- The submitted manuscript states 100 simulation replicates and 100 random train-test splits.
- `run_meta.csv` records 50 replicates and 10 site-prediction splits.
- Resolve this by either rerunning formal simulations at the manuscript-stated settings or revising the text to match the actual run metadata.

Do not overwrite submitted results unless intentionally regenerating the manuscript.
