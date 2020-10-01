<!DOCTYPE html>
<html>
	<head>
		<meta charset="utf-8">
		<meta name="viewport" content="width=device-width, initial-scale=1.0">
		<title>README.md</title>
		<link rel="stylesheet" href="https://stackedit.io/style.css" />
	</head>
	<body class="stackedit">
		<div class="stackedit__html">
			<h1 id="crystIT">crystIT - Quick Start Guide</h1>
				<h5 id="authors-clemens-kaussler-gregor-kieslich">Authors: Clemens Kaußler, Gregor Kieslich*</h5>
				<p>* <a href="mailto:gregor.kieslich@tum.de">gregor.kieslich@tum.de</a></p>
				<p>This is the README file to the python-based crystIT script, which calculates information theoretical complexity parameters as proposed by S. Krivovichev (2014) and extended by W. Hornfeck (2020). Modifications for partially occupied crystallographic orbits are included as well. It provides an accessible user interface, requiring no programming experience.</p>
			<h1 id="usage">Usage</h1>
				<p>crystIT is written in Python and uses standardized crystallographic information files (CIFs) as input. In the following, the script's package dependencies, the operation of the script and the output modes are explained.</p>
				<h2 id="package-dependencies">Package Dependencies</h2>
					<p>In addition to standard libraries such as numpy, crystIT was developed and tested in the following Python environment:</p>
					<ul>
						<li> Python 3.8.3 (available at <a href="http://www.python.org">http://www.python.org</a>) </li>
						<li> ASE 3.19.1 (Atomic Simulation Environment, more information at <a href="https://wiki.fysik.dtu.dk/ase/">https://wiki.fysik.dtu.dk/ase/</a>)</li>
						<li> Spglib 1.15.0 (more information at <a href="https://spglib.github.io/spglib/">https://spglib.github.io/spglib/</a>)</li>
						<li> PyXtal 0.0.7 (more information at <a href="https://github.com/qzhu2017/PyXtal">https://github.com/qzhu2017/PyXtal</a>)</li>
					</ul>
				<h2 id="starting-the-script">Starting the Script</h2>
					<p>Open the command window of your computer and navigate to the directory containing <code>crystIT.py</code>. Write in command line:</p>
					<pre><code>$ python crystIT.py</code></pre>
					<p>Successful startup is confirmed by crystIT's welcome message:</p>
					<pre><code>Welcome to crystIT -- A Crystal Structure Complexity Analyzer Based on Information Theory
Version 0.1, release date: 2020-09-22
Written by Clemens Kaußler and Gregor Kieslich (Technical University of Munich)
Please cite the following paper if crystIT is utilized in your work:
		Kaußler, Kieslich (2020): unpublished

Input path of .cif file or directory for complexity analysis. 's' for settings. 'e' to exit.</code></pre>
				<h2 id="output">Output</h2>
					<p>There are two modes of operation: Either, CIFs can be processed one by one in <i>single file mode</i>, or directories - possibly containing multiple CIFs - may be passed to the script in <i>batch mode</i>.</p>
					<h3 id="single-file-mode">Single File Mode</h3>
						<p>In <i>single file mode</i>, the path to a CIF is simply typed into the bash and confirmed with enter. All results are displayed in the bash after calculation, whereby the complexity nomenclature introduced by Hornfeck (2020) is applied. A sample output for K<sub>3</sub>C<sub>60</sub> is presented here:</p>
						<pre><code>------------ C:\K3C60.cif ------------
assumed formula  C20K
assumed SG       Fm-3m (225)
SG from CIF      F m -3 m (225)
lattice [A]      a: 14.24, b: 14.24, c: 14.24
angles [°]       b,c: 90.00, a,c: 90.00, a,b: 90.00
---
252.000000       atoms / unit cell
63.000000        atoms / reduced unit cell
123.000000       positions / reduced unit cell
8.000000         unique species
5.000000         coordinational degrees of freedom
--- combinatorial (extended Krivovichev) ---
2.648242         I_comb                  [bit / position]
3.000000         I_comb_max              [bit / position]
0.882747         I_comb_norm             [-]
325.733784       I_comb_tot              [bit / reduced unit cell]
0.451225         I_comb_dens             [bit / A^3]
--- coordinational (Hornfeck) ---
0.970951         I_coor                  [bit / freedom]
2.321928         I_coor_max              [bit / freedom]
0.418166         I_coor_norm             [-]
4.854753         I_coor_tot              [bit / reduced unit cell]
0.006725         I_coor_dens             [bit / A^3]
--- configurational (extended Hornfeck) ---
2.820700         I_conf                  [bit / (position + freedom)]
3.700440         I_conf_max              [bit / (position + freedom)]
0.762261         I_conf_norm             [-]
361.049612       I_conf_tot              [bit / reduced unit cell]
0.500146         I_conf_dens             [bit / A^3]</code></pre>
					<h3 id="batch-mode">Batch Mode</h3>
						<p>In <i>batch mode</i>, the path of a CIF-containing directory is typed into the bash and confirmed with enter. The results as well as warnings and error messages are compiled into a character-separated values (.csv) file which is saved as <code>batch_TIMESTAMP.csv</code> into the processed directory. Attention! With default settings, only CIFs directly present in the folder passed to crystIT are considered, subfolders are ignored.</p>
				<h2 id="settings">Settings</h2>
					<p>The settings menu is accessed by typing <code>s</code> and hitting enter.</p>
					<pre><code>Input float as symmetry tolerance 0 < x < 1      (currently 0.005).
Input int as maximum number of threads           (currently 12)
'd' to toggle between decimal separators         (currently '.').
'o' to toggle occupancy editing options          (currently False).
'r' to toggle recursive subdir scan              (currently False).
's' to toggle entropy calculation                (currently False).
'e' exit to main menu:</code></pre>
					<ul>
						<li>Input of a decimal number between zero and one changes <i>symprec</i> which defines the tolerance in cartesian coordinates for Spglib to find symmetry and simultaneously is the threshold cartesian coordinate value for identification of duplicate atom entries in the CIF: <code>|x′ − x| &lt; symprec</code>. Always use <code>.</code> as decimal separator to change <i>symprec</i>!</li>
						<li>The maximum number of threads for multiprocessing in batch mode is automatically set to the maximum number of available threads but can be adjusted by integer input.</li>
						<li><code>d</code> toggles the decimal separator between dot and comma, especially useful for German Excel users.</li>
						<li>The occupancy options, accessible by typing <code>o</code>, allow for on-the-fly occupancy editing in single file processing.</li>
						<li>By activating the recursive subdirectory scan with <code>r</code>, subfolders are scanned in batch mode.</li>
						<li><code>s</code> toggles the calculation of entropy values from information content values, according to Krivovichev (2016).</li>
						<li>Finally, the settings menu is exited with <code>e</code>.</li>
					</ul>
	</body>
</html>
