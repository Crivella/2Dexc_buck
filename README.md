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
  4. Edit the ENVIRONMENTAL_VARIABLE file by setting 
   * BIN_DIR as the path to a working version of QE (If QE is installation is set into PATH bin dir can be left empty and the script will get the correct directory using the command 'which')
   * PSEUDO_DIR as the path containing the pseudopotentials
   * TMP_DIR as the path of the temporary folders produced by QE
   * RUN_COMMAND as the mpi command for parallel execution
  5. Edit the syste.sh file to reflect your physical system and the needs of your calculation as explained by the comments in that file
  6. Run the calculation using the 'main.sh' script.  
  ```
  ./main.sh
  ```
  The script will identify already performed calculation and read their previous output instead of re-doing them.  
  If you want to redo a calculation, its previous output should be deleted.  

## PPtools
GNUplot post-processing tool for the final {PREFIX}_SAVE.dat file.
Launch these scripts from the folder containing the 'system.sh' and '{PREFIX}_SAVE.dat' files
  1. min_cd.gnu:  Plot the total energy VS buckling for every celldim specified in the LIST in 'system.sh'. Fit the data using a 'parabolic/quadratic fit' and take the total energy at the extrema for every celldim. Plot the "Total energ VS celldim(minimized)" plus its fit.
  2. min_bk.gnu :  Plot the total energy VS celldim for every buckling specified in the LIST in 'system.sh'. Fit the data using a 'parabolic/quadratic fit' and take the total energy at the extrema for every buckling. Plot the "Total energ VS buckling(minimized)" plus its fit.
  3. samedist.gnu:  Plot the "Total energy VS buckling" where the distance between atom is kept constant as the reference equilibrium distance. Also plot the fit of the "energy at the valence and conduction band extrema VS celldim"
  4. inplane.gnu:  Plot the "Total energy vs celldim" for buckling = 0 and fit it.




