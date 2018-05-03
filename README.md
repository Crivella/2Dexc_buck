# Two-dimensional exciton calculations
This collection of scripts uses Quantum Espresso(QE) to calculate the optical exciton
in two-dimensional systems, with the possibility to introduce strain and buckling.

## Requirements
  1. Working version of QE(>=6.0)
  2. Gnuplot (>=5.0)
  3. QEPPlib + c_scripts (provided in the tools directory)

## Usage
  1. Make sure to have an updated version of gcc and mpicc
  2. Compile the provided version of QEPPlib and the 4 .c programs
      ```
      cd tools
      make
      ```
      The provided makescript will perform compilation and installation in the tools/bin directory
  3. Copy the templ_* files and remove templ_ from the name
  4. Edit the ENVIRONMENTAL_VARIABLE file setting 
   * BIN_DIR as the path to a working version of QE (If QE is installation is set into PATH bin dir can be left empty and the script will get the correct directory using the command 'which')
   * PSEUDO_DIR as the path containing the pseudopotentials
   * TMP_DIR as the path of the temporary folders produced by QE
   * RUN_COMMAND as the mpi command for parallel execution
  5. Edit the syste.sh file to reflect your physical system and the needs of your calculation as explained by the comments in that file
  6. Run the calculation using the 'main.sh' script
    ```
    ./main.sh
    ```
  The script will identify already performed calculation and read their previous output instead of re-doing them. If you want to redo a calculation, its previous output should be deleted
