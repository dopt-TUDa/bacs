# BACS

BACS is a **B**ranch **A**nd **C**ut solver for **S**table Sets. It is based on the branch and cut solver [SCIP](https://scipopt.org/) and solves weighted maximal stable set problems given in DIMACS format.

## Usage

To run BACS, a SCIP installation is needed. Set the shell variable `$SCIP_PATH` to point to the SCIP main directory. Compile SCIP and BACS by using the Makefile (cmake is not yet working): `make OPT=opt` or `make OPT=dbg` to compile in optimized or debug mode. Run a single instance with `bin/bacs <file.dimacs>` to solve an instance saved in `<file.dimacs>`.

Possible command line options:

- `-s settings/<settingsfile.set>` to read settings parameters.
- `-t <time>` to set a specific timelimit.
- `-n <nodes>` to set a specific nodelimit.
- `-d <freq>` to set a specific display frequency.
- `-v <verblevel>` to set a specific verblevel.
- `-l <file.sol>` to start with an initial solution stored in solution file.
- `-m <memory>` to set a specific memorylimit.
- `-c <cutoff>` to set a specific objective cutoffbound.

Run `make OPT=opt TEST=<testset> SETTINGS=<settingsfile> test` to run a whole testset using custom settings. For this, `check/testsets/testset.test` needs to contain the paths to each instance and `check/testsets/testset.solu` can contain solution values. Predefined testsets are `short` and `minimal`.

## Repository Structure

- **check**: testsets, evaluation of (cluster) results
- **doc**: html and css for creation of documentation. Run `make doc` to create documentation.
- **graphs**: small test instances.
- **settings**: some settings files.
- **src**: source folder.

## DIMACS format

Each line starts with special symbol to define line usage (node indices range from 1 to n)

- c: comment lines (ignored)
- p: problem definition, line should start with `p edge <n> <m>` with `n` giving the number of nodes and `m` giving the number of arcs
- n: node description, followed by index and weight (`n <index> <weight>`)
- e: edge desription, followed by source, target (`e <source> <target>`)

Every line is ended with an newline character (don't forget the last line).

If at least one weight is given, all other nodes have weight 0. If no weights are given, the unweighted problem is solved (every node has weight 1).

An example file looks like:

```txt
c small example
c this could be your comment!
c
p edge 4 3
e 1 2
e 2 3
e 2 4
n 1 1
n 2 4
n 3 1
n 4 1
```

## License

Apache License, Version 2.0

## Authors

The Discrete Optimization Group, TU Darmstadt.
 - Jonas Alker
 - Annika Jäger
 - Erik Jansen
 - Marc E. Pfetsch
