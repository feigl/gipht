#!/bin/csh
cd ../gipht28.ENV_T115
sbatch -J ENV_T115 -p geoscience run_gipht_n_pairs_aci.csh
cd ../gipht28.ENV_T344
sbatch -J ENV_T344 -p geoscience run_gipht_n_pairs_aci.csh
cd ../gipht28.ERS_T072
sbatch -J ERS_T072 -p geoscience run_gipht_n_pairs_aci.csh
#cd ../gipht28.ERS_T115
#sbatch -J ERS_T115 -p geoscience run_gipht_n_pairs_aci.csh
cd ../gipht28.ERS_T179
sbatch -J ERS_T179 -p geoscience run_gipht_n_pairs_aci.csh
cd ../gipht28.ERS_T344
sbatch -J ERS_T344 -p geoscience run_gipht_n_pairs_aci.csh
cd ../gipht28.ERS_T451
sbatch -J ERS_T451 -p geoscience run_gipht_n_pairs_aci.csh


