synopsis:

		This program is for finding metal and base pair interactions.


install:
	Download at least all the files of sys and bin directory.
	Transfer the 'bin/metallic' executable file to some path folder
	and transfer all the files of 'sys' folder to some suitable location.
	We prefer the following.

		sudo cp bin/metbp.linux   /usr/local/bin/

		sudo cp sys/*   /usr/local/bin/


setup:
	For setting the environment, do the following from your bash shell
	of write the same in your .bashrc file. We prefer the following.
	Note: System will read the system files from this folder.

		export NUCLEIC_ACID_DIR=/usr/local/bin


run:

	metbp.linux [OPTIONS] [mmCIF/PDB]

example:

		metallic 1ehz.cif

bug-report:

		Parthajit Roy, 
			roy.parthajit@gmail.com

		Dhananjay Bhattacharyya,
			dhananjay.bhattacharyya@saha.ac.in 

