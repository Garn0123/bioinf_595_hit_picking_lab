# Lab 11 - Assessing Molecular Docking Results
In class this week we discussed the different components of large-scale docking,
which allow us to screen extremely large libraries of molecules (on the scale
of billions) in a relatively short amount of time with only a little bit of 
clever library enumeration and sampling tricks.

Rather than do an entire screen here, we are instead going to be looking
at molecules from a small screen I've already performed for you. The focus here 
will be on looking at molecules and making decisions from the results.

At the end of this, you'll have hand-selected a handful of molecules that you
think will make good binding interactions with our target, defended them, and 
run a Boltz-2 prediction on them to see what their predicted affinities are.

## Things to install

For this lab, you'll need to set up:
 [**ChimeraX**](https://www.cgl.ucsf.edu/chimerax/)
 
 An environment for [Boltz-2](https://github.com/jwohlwend/boltz) (though you'll just need it at the end)

## The Screen (that I ran)

### Quick note on system setup
The major pieces necessary to perform a screen with DOCK6 are:
- The protein (and optionally the cognate ligand, both usually a xtal structure)
- A collection of molecules in .mol2 or .db2 format
- An energy grid
- Docking spheres
- A DOCK6 input file

The list itself is short, but there are a lot of considerations when preparing
the system. A brief and non-exhaustive list, just for you to keep in your head 
if you're going to run your own screens later:
- Crystal structure resolution
- Crystallization conditions (and what that means for the model)
- The force field and charge model you're going to use
- Protein environment (is it a membrane protein, for example)
- Protonation states of receptor residues (especially around metal centers)
  - Also if charge distribution around metal centers will affect things
- Formal charge of the ligand
- Handling missing loops
- What spheres to use (ligand vs. pocket, and which pocket?)

Some things are straightforward, some things are kind of already solved, and what
those things are changes for every target! Keeping good notes and making sound
decisions you're comfortable defending are usually the name of the game.


### How the molecules were docked for this lab
For this lab, we're going to screen against [5ZV2](https://www.rcsb.org/structure/5ZV2),
which is FGFR-1 in complex with ligand lenvatinib. This is a cancer target bound
to an anti-cancer therapeutic. 

The system was prepped using a similar methodology to the Rizzo Lab's virtual
screening protocols (https://github.com/rizzolab/DOCK6_Screening_Protocols).

The library used for this screen was a random collection of 500k drug-like
molecules. All molecules were docked using the FLX protocol of DOCK6 (not the
HDB method discussed in class) across 240 cores on an HPC. The total time to
screen was ~a day.

I then took the top 5000 molecules by Grid Score from that screen (so the top 1%)
and minimized them in Cartesian space. We used grids to save a bunch of time, 
but now that we don't have to sample a bunch of conformations we'd prefer the
more accurate pairwise Cartesian pose over the grid-specific pose to make it
more comparable with the Xtal ligand.

I then used RDKit (via DOCK6) to calculate molecular features and then rescored 
mols using the following scoring functions:
- Continuous Score (pairwise LJ + Coulombic) (negative better)
- Footprint Similarity (FPS) Score (Measures how well energy profile of ligand is to a reference, ie. how different are the interactions per residue of the protein) (Higher = less overlap)
- Hungarian Matching Similarity (HMS) Score (RMSD-like term that can compar dissimilar molecules as well as account for ring flips) (-5 is perfect overlap, infinite is no overlap)
- Volume Overlap Score (VOS) (How well the ligand fits into the reference ligand volume)

Remember that rescoring doesn't change the pose! These are the 'best' energetic
poses, but we just calculated scores for that pose after the fact.

The similarity-based scoring functions compare features of the molecule against 
a reference - in this case, against the crystal lenvatinib. Sometimes, if we 
have a molecule we know already binds, we might be interested in matching similar
energetic features. 

Things like FPS score can be pretty useful to explore the pocket landscape, 
but I am not providing the individual FPS files here to actually construct a
pairwise visualization. Just take them as a matching measure rather than an 
individual measure and we can talk more about it if you wish.

For more information about the individual scoring functions than my notes above,
please reference the [DOCK6 user's manual](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm)


So quick recap:

500k mols docked -> Top 5k by grid score -> Minimized top 5k -> Rescored top 5k


## The Lab
### Selecting 10 molecules in Chimera
Open the docked molecules (~/docked_molecules/5ZV2_5k_cartmin_rescored.mol2) in
a text editor. You'll see a series of MOL2 blocks broken up by a "header" that
describes each molecule. Each header entry is prepended by a bunch of pound signs,
making them easy to extract via a Python script or other such tool. Everything 
with "RD" in it is an RDKit calculated feature, and anything with "desc" in front
of it is a scoring function related value. 

This contains all the information from the docking steps in the previous section.
If you open this file in ChimeraX, it'll open a ViewDock window that will let
you page through the molecules one-by-one and see their scores side-by-side.


Open up the following files in ChimeraX:

- `~/system_files/5ZV2.rec.clean.mol2`
- `~/system_files/5ZV2.lig.am1bcc.scored.mol2`

If you're new to Chimera, take this opportunity to explore it a little bit. Before 
what you need to turn in at the bottom of this lab I'll put some tips and tricks
that can help you get started.

I would calculate H-Bonds between your protein and the xtal ligand, just to get
a baseline for what we're looking at. You can do this via:

Tools > Structure Analysis > H-Bonds

I usually select my ligand or protein and calculate intermodel H-Bonds and 
reveal the atoms of H-bonding residues.

Take a look at the score of the molecule and how it packs into the site (you
may also want to look at the protein surface).

Now open the docked molecules and calculate H-Bonds between all of the ligands
and the protein. Page through the molecules and take a look at their H-Bonds,
packing, and their scores to get a feel for the structures and potential interactions.

There are a ton of molecules, and not all of the numbers are useful! This is a 
situation we often find ourselves in when we're thinking about these systems.
Think of different ways to sort the data that might give you molecules from 
different classes. Don't just trust the score! 

I also don't necessarily expect you to look at all 5000 molecules, but you should 
get a fair feel for what's going on and think through how you would work through 
the list.

From what you know about binding interactions and what these numbers mean, pare
down this list to your top 10 molecules that you would suggest to carry forward
as if you were collaborating with an experimentalist. ***Write a note for each molecule
that says why you chose it for the top of your list. Include how you got there -
how did you sort the list, what specific interactions, etc.*** 

Extract those 10 molecules from the file (which you can do from ViewDock) and 
take note of the SMILES strings for the Boltz-2 step. Rank order your selections
by whatever criteria you wish, just put them in an order of "best binder" to 
"worst binder" of your set.

### Co-fold your molecules with Boltz-2
In the `~/boltz_2` directory you'll find the yaml that you'll need to assess
the 5ZV2 target with Boltz-2. Run Boltz with that yaml file. 

Using that yaml as a base, do the same for each of your 10 molecules.

Take a look at how different the pose is in Boltz-2 (preferably with the
`matchmaker` command in Chimera). Note down the predicted affinities for your
molecules.

## Chimera Starter Tips
You can use Ctrl + Click to select things in the window, and things you click
on in the GUI will affect what you have clicked. You can clear your selection
by Ctrl + Clicking the background or in the menus go to Select > Clear.

If you select an atom you can press up or down and it'll move up to the next
largest selection (from atom -> substructure -> entire model).

You can change the color of your selection with the Actions > Color window.

You can hide non-polar hydrogens by doing: 
`Select > Chemistry > IDATM Type > HC` and then `Actions > Atoms/Bonds > Hide`

Most everything visualization-wise can be accessed from the top toolbar - don't
be afraid to press buttons!

You CAN undo a lot of commands/actions with Ctrl + Z.



## The Questions
1. What was your process for trying to select molecules? Did you focus on any
specific values in your search? 
2. This was a significant amount of data, even for "only" 5000 molecules. At scale,
even 0.1% of the library could be several orders of magnitude greater than that 5000.
What might be some ways you can think of to make that kind of volume of data 
tractable to assess?
3. Compare your own ordering of your molecules to the affinities from Boltz-2. 
Did those two lists match orders? Give some thoughts on why or why not.
4. You probably noticed differences when comparing your Boltz-2 predictions to the xtal
pose, particularly in placement of some protein residues. DOCK6 is a rigid docking
model - the protein doesn't move at all. Does this limitation make this a
reasonable model for this application? Why or why not? Can you think of any other
limitations to this methodology? Would addressing those limitations change the
usability of the model in a negative way?



## What to turn in

- A file containing your 10 selected molecules (selected_mols.mol2)
- Your 10 Boltz-2 models
- A report with the following information:
	- The predicted affinity values for all 10 Boltz-2 models
	- The confidences for the Boltz-2 predictions
	- Why you selected your 10 molecules prior to refolding with Boltz-2
- The answers to `The Questions` above.
