import os

DoSUBMIT = True

PNFS_BASE = '/pnfs/icarus/scratch/users/rvizarr/toy_ensemble'
gundam_Tar = f'{PNFS_BASE}/tars/local_gundam.tar'
config_Tar = f'{PNFS_BASE}/tars/local_config.tar'

NToys = 1000           # 1 → 5 → 50 → 1000
OutputDest = f'{PNFS_BASE}/output'

FitterConfigPath = 'RunConfigs/dpT/Asimov_StatSyst/config_Fitter_FakeData_dpT.yaml'

jobsub_cmd = f'''jobsub_submit \\
-G icarus \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
--lines '+FERMIHTC_AutoRelease=True' \\
--lines '+FERMIHTC_GraceMemory=1000' \\
--lines '+FERMIHTC_GraceLifetime=3600' \\
--disk 20GB \\
--memory 4GB \\
--expected-lifetime 6h \\
-N {NToys} \\
"file://$(pwd)/grid_executable.sh" \\
--dest {OutputDest} \\
--gundamTar {gundam_Tar} \\
--configTar {config_Tar} \\
--FitterConfigPath {FitterConfigPath}
'''

print(jobsub_cmd)
if DoSUBMIT:
    os.system(jobsub_cmd)
