# README improvements — pre/post publication

## Blocking

- [ ] **Glimmer training file — add a download source.**
  Line 145 currently says *"A download link will be provided in a future release"* — but
  this IS the release. `training-file_refseq.icm` is required to run the pipeline and there
  is currently no way for a user to obtain it. Either host it on Zenodo / a GitHub release
  and add the URL to `setup_databases.sh` + the README, or document clearly where it comes
  from so users can reproduce it themselves.

## Should fix

- [ ] **Inconsistency between copy instruction and run command.**
  The Configuration section tells users to `cp config.yml my_config.yml`, but every run
  command still shows `-params-file config.yml`. Either drop the copy suggestion, or update
  all run commands to use `-params-file my_config.yml`.

- [ ] **Add disk space estimate to the Database setup section.**
  PHROGs + Pfam-A + ECOD are collectively ~15–20 GB. A single sentence saves users from
  starting a multi-hour download and running out of quota on a cluster.

- [ ] **Pin HH-suite version in `environment.yml` and the Software versions table.**
  Currently listed as *"latest via bioconda"*, which is not reproducible for a published
  pipeline. Pin the version that was validated and reflect it in the README table.

## Nice to have

- [ ] **Add a Quick Start section** near the top — three commands (clone, setup DBs, run)
  for users who already know Nextflow and just want the gist.

- [ ] **Add typical runtime / memory guidance.** HHsuite searches against Pfam/ECOD can be
  memory-hungry. A note on expected wall time and memory requirements would help users size
  their HPC job requests.
