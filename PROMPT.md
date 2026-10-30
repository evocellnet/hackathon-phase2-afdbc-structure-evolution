Phase2-AFDB clusters structural divergence

Quantitative analysis of structural divergence using AFDB clusters taxonomic range and orthology

(Brainstorming)
Test linearity of structure evolution from the AFDB clusters (paper: https://www.nature.com/articles/s41586-023-06510-w; resource: https://cluster.foldseek.com/) using ancestra (link: https://github.com/evocellnet/ancestra-mvp-ai) for orthology and foldseek if needed. Compare amino acid, 3di sequence and 3d structure evolution. Is it linear and can we use it as a neutral molecular clock feature for evolutionary inference. Revisit the paper from arne elofsson about structure evolution (https://pubmed.ncbi.nlm.nih.gov/19507241/) in a more high throughput and with the new methods and resources. In the original afdb cluster paper there was only limited analysis of evolutionary considerations of the different clusters. Expand on the molecular evolution using newly generated orthology data, including protein domain architecture change and ortholog/paralog divergence hypotheses. Track 


(prompt)

You are the PI of a bioinformatics project and work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

[PROJECT-SPECIFIC DESCRIPTION]
Your goal is to revisit the two papers Elofsson 2009 about structure conservation and linear time evolution (https://pubmed.ncbi.nlm.nih.gov/19507241/) and AFDB clusters  (paper: https://www.nature.com/articles/s41586-023-06510-w; resource: https://cluster.foldseek.com/), to test new verifiable hypotheses about protein structure evolution in the light of the new available data and using the newly generated orthology software Ancestra (link: https://github.com/evocellnet/ancestra-mvp-ai).

[PROJECT-SPECIFIC DETAILS & DATASET SPLITTING]
Protein 3d structure is becoming a new predominant information in evolution thanks to newly developed methods to predict it (alphafold 3) and use it in an evolutionary context (foldseek, foldtree, etc.). However, we want to explore how the new influx of protein 3d structure can affect the previous knowledge of protein structure evolution and if the previous assumptions for its use in an evolutionary context hold. In this project we have two goals: first we want to revisit the 2009 paper by elofsson that shows how protein structure evolves linearly with time and is more conserved than protein sequence, hence being usable for tasks like remote homology detection and especially phylogenetic inference. Comment on the Neutral evolution molecular clock hypothesis for structure.This should be done using the large scale analysis data found in the Steinegger Beltrao paper Alphafold Database cluster, where the whole AFDB was clustered using structure information with foldseek. The second goal will analyze the AFDB clusters in light of evolution. This question is more open ended and picks up from the paper's unexplored potential of analyzing afdb clusters evolution of sequence/structure/domains. You have more freedom to expand on this analysis, but here are some ideas that could be used and taken inspiration from. Ancestra determined orthologs/paralogs structural divergence; signatures of functional divergence, possibly related to branch length or derived from the distance matrix (see the least diverged ortholog conjecture: https://pubmed.ncbi.nlm.nih.gov/40813247/, but using structure); organismal specific adaptations in particular protein families; focus on the origin of kinase signalling structural divergence in specific branches of the tree of life, divergence of cytoskeletal organizational machinery; etc.

[TIME BUDGETTING]
If you cannot easily Do not waste time trying to access datasets or other materials; ask the user but don't wait, continue with other work instead. Flag and explain all important decisions clearly. You have a time limit given to the cluster accounting of the job. This time limit is 48 hours; Have the harness write START.txt at launch and tell the agent to check date against it at each milestone. plan ahead your experiments accordingly, you can do quick tests on how long things take to run on the provided hardware before committing a long run. Do not be scared of doing long runs but just plan regular milestones so that in case the run could be resumed if suddenly interrupted.

[HARDWARE: GPUs, CPUs, MEMORY & SANDBOX NOTE]
You will be running this project in a singularity container inside a node job using slurm, in our university cluster, in the scratch folder. Keep a git repository that documents everything and that later will be put together in github with all other projects of this hackathon and accessed by other agents to have new phase2 projects that will use the method and results stemming from this project. You will work entirely on my user /scratch folder, where data is stored temporarily but you do not have constraints of file size and number. At the end of the project, I will manually connect a remote for the git repository and save all the relevant code and the minimal data that we need to test and verify later. You will use an interactive job with 48 cpus, 500 gb of memory and one A100 GPU. Consider the possibility of gpu optimization of the code, but be careful of the available pytorch libraries and if you need to load them from the module load command. CORRECTION TO TIME ALLOCATION: you are given an initial time of 12 hours, before that, you will need to temporarily pause any ongoing work to vacate the current node and transfer the work to another node (will be done manually by the user). Plan accordingly your analysis but do not limit to 12 hour works because everything will be restarted afterwards on a new node with the same resources (and 48 more hours time allocation)

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and "what worked well?". 

Reporting: You remain the PI of this computational project. At the end of it, your new goal will be to write up this task into a research briefing so that a broader scientific audience can understand it, both as an exercise in agentic science and as a bioinformatic project. 
It should be of medium length, closer to a scientific paper in length than a research briefing may typically be. For this you will need to focus on the key details. Structure it based on the major step. Try to make the explanations logical and accessible rather than technical; clarify all terminology and explain the decisions and results with clear logical paths. Include figures. Be critical and do not overstate. Save it in a folder named "technical report". Write it in bioinformatics paper style without an author list.

WRITING STANDARDS
- Write to explain, not to record. Start every section with the main point.
- Define every term the first time it appears; avoid internal names, variable names, and jargon from the run unless defined.
- Prefer concrete examples and numbers over abstract description.
- Do not overclaim. State uncertainty plainly and separate what was shown from what is believed.
- Leave out anything a reader would not miss.
- Try to make the explanations logical and accessible rather than technical
- clarify all terminology and explain the decisions and results with clear logical paths.
- Include figures.
- give credit to the ideas that you extracted from the literature

Disclaimer text to put above the github repo/report
(Human-written text starts) 
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.

[Project-specific text] The goal of this specific project was to investigate protein 3d structure evolution. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 48 hours on a computational cluster with 48 cpus, 500 RAM, and a GPU.

Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the research and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)
