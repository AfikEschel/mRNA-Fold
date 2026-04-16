mRNA secondary structure prediction using
utility-scale quantum computers
Dimitris Alevras∗, Mihir Metkar‡, Takahiro Yamamoto†, Vaibhaw Kumar∗, Triet Friedhoff∗,
Jae-Eun Park∗, Mitsuharu Takeori†, Mariana LaDue∗, Wade Davis‡, and Alexey Galda‡
∗IBM Quantum, New York, USA
†IBM Quantum, Tokyo, Japan
‡Moderna, Cambridge, USA
Abstract—Recent advancements in quantum computing have genetics, aiding in the design of mRNA-based drugs, and
opened new avenues for tackling long-standing complex combi- developing gene therapies.
natorial optimization problems that are intractable for classical
The problem of RNA secondary structure prediction with
computers.PredictingsecondarystructureofmRNAisonesuch
pseudoknots is known to be NP-complete [5]. State-of-the-
notoriously difficult problem that can benefit from the ever-
increasingmaturityofquantumcomputingtechnology.Accurate art classical approaches for the secondary structure prediction
prediction of mRNA secondary structure is critical in designing of mRNA have primarily relied on dynamic programming
RNA-basedtherapeuticsasitdictatesvariousstepsofanmRNA algorithms, such as the Zuker algorithm [6] and its widely
life cycle, including transcription, translation, and decay. The
used implementations, MFold [7] and ViennaRNA [8], to
current generation of quantum computers have reached utility-
name a few. These methods calculate the minimum free
scale, allowing us to explore relatively large problem sizes. In
thispaper,weexaminethefeasibilityofsolvingmRNAsecondary energy (MFE) structures by systematically evaluating all pos-
structures on a quantum computer with sequence length up to sible base pairings and applying thermodynamic models to
60 nucleotides representing problems in the qubit range of 10 estimate their stability [9]. Subsequent enhancements have
to 80. We use Conditional Value at Risk (CVaR)-based VQE
incorporated stochastic context-free grammars (SCFGs), as
algorithm to solve the optimization problems, originating from
seenintoolslikeInfernal[10],whichalignsequencestoRNA
the mRNA structure prediction problem, on the IBM Eagle
and Heron quantum processors. To our encouragement, even structureprofilestopredictconservedstructuresacrossrelated
with “minimal” error mitigation and fixed-depth circuits, our sequences.Theintegrationofmachinelearningtechniqueswas
hardwarerunsyieldaccuratepredictionsofminimumfreeenergy implemented in algorithms like ContextFold [11] that utilize
(MFE) structures that match the results of the classical solver bothlocalandglobalcontextualfeaturestoimproveprediction
CPLEX. Our results provide sufficient evidence for the viability
accuracy. Despite these advancements, classical approaches
of solving mRNA structure prediction problems on a quantum
computer and motivate continued research in this direction. arestillconstrainedbythecombinatorialexplosionofpossible
structures, lack of generalization to novel RNAs [12], and the
approximations required to make the problem tractable [13].
I. INTRODUCTION
As a result, there is a growing interest in leveraging quantum
Messenger RNA (mRNA) plays a pivotal role in the cen- computing to overcome these limitations and achieve more
tral dogma of molecular biology, acting as the intermediary accurate predictions of mRNA secondary structures.
between the genetic blueprint in DNA and the functional ClassicalhardnessofmRNAsecondarystructureprediction
molecule, protein. The functionality of mRNA is inherently naturally positions it as a compelling application for quantum
linkedtoitssecondarystructure,whicharisesfromintramolec- computers. Despite still being relatively scarce, recent devel-
ularbasepairingandfoldingpatterns[1].Understandingthese opments have shown promising advancements in this area. In
structures is crucial for elucidating the mechanisms of gene Refs. [14], [15], quantum annealing-based hardware was used
regulation, translation, and degradation kinetics, and inter- tominimizetheenergyfunctionthatcorrespondstothefolding
biomolecular interactions [2]. Additionally, the emergence energy landscape of mRNA molecules. However, quantum
of mRNA-based therapeutics has ignited a paradigm shift annealers are not universal quantum processors, which limits
in the landscape of medicine, offering novel solutions for the set of algorithms amenable to this technology compared
the treatment of various diseases [3]. This transition from to the universal gate-based quantum computers. In another
understanding mRNA’s fundamental biological role to har- recent paper [16], the Quantum Approximate Optimization
nessing its therapeutic potential marks a significant milestone Algorithm (QAOA) [17] was used to solve the secondary
in biomedical research and holds immense promise for rev- structure prediction problem for toy-model RNA sequences
olutionizing healthcare practices. However, predicting mRNA using 12 qubits on a simulator and 4 qubits on a gate-based
secondary structures is a formidable combinatorial problem quantum processor.
due to the vast number of possible configurations and their AlgorithmssuchasQAOAandVariationalQuantumEigen-
associated energy states [4]. Solving this problem accurately solver (VQE) fall under the class of variational quantum
and efficiently is essential for advancing our knowledge in algorithms that can be implemented on the current generation
4202
yaM
03
]hp-tnauq[
1v82302.5042:viXra
of gate-based quantum computers. Guided by a classical results and an outlook.
optimization feedback loop, the goal of such algorithms is
to find the quantum state that exhibits sufficient overlap with
the ground state of the optimization problem. Success of
such schemes depends heavily on the choice of a suitable
ansatz, objective function, and efficient classical optimizer, as
highlighted in the literature. However, most of such studies
areeitherbasedonnumericalsimulationsorexaminerelatively
smallproblemsonquantumhardware,withonlyafewnotable
exceptions [18], [19], [20], [21], [22].
As universal quantum processors surpassed one hundred
qubits, we have entered the “utility” era of quantum com-
puting [23], where quantum computers can perform reliable
computations on problem sizes beyond the range that can be
handledbybrute-forceclassicalschemes.Encouragedbythese
developments,inthisworkweexaminethefeasibilityofsolv-
ing relatively longer mRNA sequences on a universal utility-
scale quantum processor. The goal here is to demonstrate that
variational quantum algorithms can be fully and successfully
executed on such hardware, with hundreds of quantum circuit
calls during the hyperparameter training step along with basic
error mitigation, yielding accurate results that match those
obtained with classical optimization schemes.
In this work, we use the Conditional Value at Risk Vari-
ational Quantum Eigensolver (CVaR-VQE), proposed as a
modification to the traditional VQE algorithm, that employs
CVaR values of the energy eigenfunctions as the objective
function, to achieve better convergence [24]. Additionally,
we choose Nakanishi-Fujii-Todo (NFT) [25] as a classical
optimizeralongwiththehardware-efficient“two-local”ansatz.
The combination of the NFT optimizer and two-local ansatz Fig. 1: An RNA secondary structure illustrating five kinds of
has been shown to be relatively robust against hardware noise secondarystructuralelements:hairpinloops(blue),bulge(yel-
and exhibits better convergence properties [26]. low), helix (black), internal loop (green) and multi-branched
First, we formulate the secondary structure prediction loop (cyan) and a pseudoknot (pink) [28].
problem as a binary optimization problem. We then use
CPLEX [27] to solve the formulated problem. We provide
numericalevidencethatastheproblemsizeincreases,classical
II. PROBLEMDEFINITION
schemes like CPLEX show an extremely unfavorable runtime RNA commonly is a polymer of four nucleotides (bases):
scaling.Next,weruntheseproblemsusingbothquantumhard- adenine (A), cytosine (C), guanine(G), and uracil (U). Un-
ware and simulations, and present results on RNA sequence der physiological conditions, RNA adopts compact structures
lengthsrangingfrom15to42basesrepresentingoptimization driven by hydrophobic, base-stacking, and base-pairing in-
problems in the range of 10 to 80 qubits. To the best of our teractions [1]. Typically these interactions facilitate the for-
knowledge, our study is the first to examine RNA sequences mation of an A-form helix containing canonical Watson-
up to 42 nucleotides on a quantum processor. Crick-Franklin base-pairs. The folded RNA conformation can
The paper is organized as follows. In Section II, we offer a be depicted as secondary structures, comprising various sec-
formaldefinitionoftheproblemandgiveabriefdescriptionof ondary structural elements that arise from base-pairing and
thestructuresthatprovidemoreorlessstability.InSectionIII, base-stacking interactions (see Fig. 1). Understanding RNA
we briefly describe the current classical methods and describe secondary structure is crucial for gaining insights into its
state-of-the-art classical algorithms, with focus on dynamic functionanddesigningRNA-basedtherapeutics[2].TheRNA
and mathematical programming approaches. We also provide secondary structure prediction problem is the problem of
our formulation of the problem as a mathematical program- finding the most stable folding of the sequence of bases [30].
mingproblem.SectionIVfocusesonthequantumformulation The exact characteristics of this “optimal” folding are not
of the problem and describes the CVaR-VQE method used in completely understood, making approximate models well-
this work. Section V contains the results of our experiments, suited for this purpose. In fact, in certain cases, there may
using both classical simulators and quantum processors. Fi- exist multiple optimal structures, underscoring the value of
nally, in Section VI, we conclude with a summary of our methods that can produce numerous high-quality solutions.
• AQuadraticUnconstrainedBinaryOptimization(QUBO)
problem can be converted to the problem of finding the
lowestenergyeigenstateofanIsingHamiltonianandthus
can be solved by a quantum algorithm;
• The secondary structure prediction problem can be for-
mulated as a QUBO.
Previouseffortshavefocusedondevelopinganobjectivefunc-
tion to minimize over binary variables, with subsequent opti-
mization carried out on a quantum annealing hardware [14],
[15]. In this work, we use mathematical programming to
represent the problem in binary variables and then convert the
resulting program to QUBO, which we solve on a universal
quantum computer.
A. Dynamic programming algorithms
One approach to approximate the “optimal” folding is to
maximize the number of paired bases. This can be accom-
Fig. 2: Example of an RNA sequence and its predicted plished by the dynamic programming algorithm of Nussi-
MFE structure demonstrating its folded state. Optimal folding nov[30].Thealgorithmisbasedonthefactthatthemaximum
calculated by ViennaRNA RNAfold server. Graphic generated numberofpairedbasesforsequencei,...,j canbecalculated
with ViennaRNA forna server [29]. from the maximum number of pairs for smaller length se-
quences.Thereareonlyfourpossiblewaystoconstructnested
base-pairs from smaller sub-sequences [37]. It is therefore
Accordingtothermodynamicprinciples,formationofstem- possible to recursively calculate the maximum number of
loop structure is feasible if the entropic costs of closing loops pairs. Even though maximizing the number of paired bases
and ordering nucleotides in a helix are offset by free energy is not sufficient, Nussinov’s algorithm provides the basic
releasedduetobase-stackingandpairinginteractions[31].As mechanism of subsequent dynamic programming algorithms,
aresult,longerhelices,orstems,aretypicallymorestabilizing suchasZuker’salgorithm[38];onehastoreplacethenumber
for secondary structures, whereas loops and single-stranded of maximum pairs in Nussinov’s algorithm with a scoring
regions tend to introduce instability. To predict the most function (such as MFE) and calculate the score for each sub-
favorable folding, a computational model must determine the sequence. The best approximation to what “optimal” means
structure with the lowest free energy (MFE) compared to the is the one that results in the best score (e.g., MFE) for the
unfolded state, which is the conformation that RNA is most sequence.
likely to adopt under physiological conditions [31]. Dynamic programming algorithms do not account for all
structures, especially a specific type of structure known as
III. CLASSICALAPPROACHES a pseudoknot, shown in Fig. 1. While not present in all
mRNA sequences, and particularly relatively short sequences
The most common way of predicting secondary structure considered in this work, pseudoknots are crucial for the
is to minimize the free energy change upon moving from an function of several important RNA elements, including reg-
unfolded state to the most thermodynamically stable folded ulation of translation and splicing [39] and the binding of
state, see Fig. 2 [32]. This secondary structure which has small molecules [40], [41], [42]. However, for many practi-
the lowest possible energy is defined as the MFE structure. cal purposes, approximation methods that do not explicitly
Because MFE reflects the most thermodynamically stable account for pseudoknots can still provide reasonable and
structure, many approaches for predicting secondary structure efficient predictions of mRNA secondary structures. The cur-
concentrate on identifying structures through MFE minimiza- rent classical algorithms utilized in software achieve 65-73%
tion. However, this problem, viewed as an optimization prob- accuracy in predicting secondary structures observed in the
lem, is generally NP-hard. There are approximation schemes laboratory [43].
designed to optimize simplified models by excluding specific
structureslikepseudoknots(seeFig.1)andfocusingonlyona B. Mathematical programming models
subsetofpossiblestructures[30],[33].Theseincludedynamic Amongclassicalschemes,mathematicalprogrammingmod-
programmingalgorithms(describedbelow),alongwithgenetic els and integer linear programming, in particular, are often
algorithms[34]andmachinelearningalgorithms[35],toname used to solve the RNA folding problem [44]. Although these
afew.SeeRef.[36]forareviewofclassicalmethodsforRNA modelsendupbeingverylargeinthenumberofvariablesand
structure prediction. constraints, the classical solvers available today are sophisti-
The potential of quantum computing to tackle this problem cated and powerful enough to handle problems of reasonable
stems from two key aspects: size. A significant benefit of integer linear programming lies
initsflexibilitytoaccommodatemodificationstotheproblem, • Each variable (quartet) has a free energy as given in
such as incorporating new structural motifs like pseudoknots, Ref. [45].
without disrupting the method, unlike dynamic programming • Consecutive quartets are preferred since they provide
algorithms. more stability. We model this, by adding the product of
Severalwayshavebeenproposedtorepresentthesecondary thecorrespondingvariablesintheobjectivefunctionwith
structure prediction problem, defining a variable to be: a reward.
• a pair between bases i and j represented as (i,j) • Certain structures are discouraged. In this work we con-
• a stack of pairs (i,j),(i+1,j −1),...,(i+k,j −k) sider the (UA) pair ending penalty as given in Ref. [45].
represented as (i,j,k) We implement this by adding variables that form such
• astackoftwoconsecutivepairs,alsoknownasaquartet, ending pair to the objective function with a penalty.
represented as (i,j,i+1,j−1). For brevity of notation, we use q to denote the ith quartet
i
Using a pair between i and j as a variable works well in variable. Let Q be the set of valid quartets, QC the set of
solvingthemaximumpairproblem(likeNussinov’salgorithm) all pairwise crossing quartets, QS(q i ) the set of quartets that
butmakesitdifficulttomodelmorecomplexrequirementsand can be stacked with quartet q i , and QUA the set of stacked
constraints.Intermsofscaling,thenumberofvariablesofthe quartets that have a (UA) end pair. The sets Q, QC, QS(q i )
problem with respect to the sequence length scales best when and QUA, are determined based on the given sequence in
using quartets and worse when using pairs. a pre-processing step. Each quartet q i ∈ Q is formed by two
We use mathematical programming to model the prob- consecutive(stacked)pairspi 1 andpi 2 .TheMFEe qi forquartet
lem as a binary optimization problem. Each of the po- q i istheMFEofp 1 followedbyp 2 [45].Twoquartetsq 1 ,q 2 ∈
sitions in an RNA sequence 1,...,n consists of one of Qthatarestackedarerewardedbyarewardr.Ifastackends
the bases U,A,C,G. A pair between two bases in po- in a (U,A) pair there is a penalty p. The objective function is
sitions i and j in the sequence is allowed if the cor- thus:
responding bases form a valid pair. The valid pairs are: (cid:88) (cid:88) (cid:88) (cid:88) (cid:88)
min e q +r q q +p q (1−q ).
{(AU),(UA),(CG),(GC),(GU),(UG)}. The variables are
qi i i j i j
defined for each quartet [44] which is two consecutive pairs
qi∈Q qi∈Qqj∈QS(qi) qi∈Qqj∈QUA
(1)
(also referred to as stacked pairs):
The constraints of the problem implement the requirement
(cid:26) 1 if (i,j), (i+1,j−1) are made that if two quartets q 1 ,q 2 ∈ Q are crossing (i.e., contain
x(i,j,i+1,j−1)=
0 otherwise crossing pairs) they cannot be both selected:
These variables are generated if the pairs are valid. The q i +q j ≤1, ∀q i ,q j ∈QC. (2)
following data pre-processing steps are performed to prepare
The resulting problem is a quadratic binary optimization
the data for the formulation:
problem. We use CPLEX as a classical solver to solve this
• check if two pairs can create a quartet (can nest) problem and create a baseline. As noted above, for larger
• check two pairs for having the same base sequencesthenumberofconstraintsinthisformulationisquite
• checks two pairs for crossing large.However,modernclassicalsolvers,suchasCPLEX,can
• create a list of quartets (variables of the problem) deal with the size of the problem. In our experiments, we
• check two quartets for crossing generated10,000sequencesoflengthvaryingbetween15and
Two quartets cannot be both selected if there are crossing 60 nucleotides, which resulted in quadratic binary problems
pairs between them. Also, each base can be paired only with withupto327variablesandupto38,946constraints.CPLEX
one other base. These conditions result in constraints of the was able to find optimal solutions to all instances within a
form: few seconds per instance. In Fig. 3 we plot the time CPLEX
takes to solve a problem, measured in internal clock ticks, as
x(i,j,i+1,j−1)+x(k,ℓ,k+1,ℓ−1)≤1.
a function of length of RNA sequences. Our dataset contains
Since the variables are binary, the constraint above is equiva- over 600 problem instances per sequence length The box
lent to requiring that: plot suggests that the time to solve the hardest instances
grows exponentially with sequence length. As a result, in
x(i,j,i+1,j−1)x(k,ℓ,k+1,ℓ−1)=0,
the mathematical programming formulation of the problem,
whichcanbeincludedintheobjectivefunctionasaproductof classical solvers will struggle to find a solution within a
the two variables with a penalty, thus preventing both of them reasonable time frame even for sequences of a few hundred
tobe1.Wenotethatdependingonthealgorithmusedtosolve nucleotides.
theproblem,wemaywanttokeepone(withconstraints)orthe
IV. QUANTUMAPPROACH
other (without constraints) formulation. In our experiments,
A. QUBO formulation
we utilize CPLEX, which demonstrates superior performance
when explicit constraints are included in the formulation. The mathematical programming formulation described
The objective function has the following components: above is a quadratic constrained binary optimization problem.
2500
2000
1500
1000
500
0
15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60
Sequence length
skcit
kcolC
|ψ(θ)⟩ is used to compute a suitable objective. The goal is to
identify the optimal parameters θ that minimize the objective
function, typically achieved through a classical optimization
scheme. We use here a CVaR-based objective that is defined
as the average of lower α-tail of the energy distribution of
the sampled bitstring. More formally, if the set of bitstring
energies sampled and sorted in increasing order is denoted by
{E } where i∈[1..K], then the CVaR(α) is defined as
i
⌈αK⌉
1 (cid:88)
E . (4)
⌈αK⌉ i
i=1
Whenα=1,CVaR(α)correspondstotheexpectationvalue
⟨ψ(θ)|H |ψ(θ)⟩ of the problem Hamiltonian H , whereas
p p
when α → 0, only the minimum among sampled energies
is considered. Employing an α value in the range (0,1] is
expected to improve the convergence behavior of the classical
Fig. 3: Time scaling (internal clock ticks) for CPLEX as optimizer [24], [48]. More recently, authors in Ref. [20]
a function of sequence length for 10,000 random RNA se- indicated CVaR to be robust against quantum hardware noise
quences. and provided bounds on the noise-free expectation values,
placed by CVaR based on the noisy samples.
C. Classical Optimizer
The corresponding problem Hamiltonian can be obtained by
transforming the constrained problem to a QUBO problem. Recently, several studies [49], [50] have examined the
Other approaches [15], [14] formulate the QUBO directly performance of classical optimizers for variational schemes,
and train the model to find the appropriate penalties and with an emphasis on robustness of such optimizers on the
rewardsfortheobjectivefunction.HereweusetheQiskit[46] noisy near-term quantum hardware. In this work, we employ
available methods to convert the problem to a QUBO and the NFT algorithm [25] as the classical optimizer which
calculate the penalties for the constraints. was found to be relatively better under hardware noise [26],
As previously mentioned, the constraints in Eq. 2 can be [51]. NFT is a gradient-free classical optimizer that falls
incorporated into the objective function (Eq. 1) through the under the class of sequential minimal optimization algorithms
introduction of a slack variable and a penalty term. However, such as Rotosolve [52], [53], [54], [55], [56]. Under this
due to the form of the constraints and the binary nature of scheme, each parameter θ i is optimized sequentially. NFT
the variables, the introduction of additional slack variables uses the fact that for parameterized circuits with parameters
is unnecessary. After adding a penalty t for the constraints, sitting only on the single-qubit Pauli rotations (and fixed
the final objective function of the resulting QUBO takes the two-qubit gates), the expectation of an observable ⟨O⟩ can
following form: be expressed as an univariate sine function of θ i with three
unknown coefficients. Under the scheme, this sine function
(cid:88)
min e qi q i is constructed by evaluating ⟨O⟩ at three different values of
qi∈Q
(cid:88) (cid:88)
θ
i
, while keeping rest of the parameters fixed. The optimal θ
i
+r q i q j corresponds to the angle that minimizes the sine function. For
q
(cid:88)
i∈Qqj∈
(cid:88)
QS(qi) convergence guarantees, unbiased estimate of ⟨O⟩ is needed.
+p q (1−q ) (3) However, numerical simulations conducted in Ref. [25] found
i j
qi∈
(cid:88)
Qqj∈QUA NFTrobustunderfinitesamplingaswell.Inthiswork,weuse
+t q q , CVaRvaluesasaproxytotheexpectationvalues.Basedonthe
i j
MPS-simulator VQE runs, we notice, NFT usually achieving
qi,qj∈QC
s.t. q ∈{0,1} ∀q ∈Q. convergence faster when CVaR with α<1 is used as a proxy
i i
to the expectation (α=1).
B. CVaR-based VQE
D. Ansatz
In this work, we use Conditional Value at Risk (CVaR)-
based VQE [24], [47] to solve the QUBO problems arising We use here the “two-local” ansatz involving single-qubit
from the formulation discussed in Section IV-A. VQE is a Pauli-Y rotations (e−iϕY) and two-qubit control-Z gates, as
hybrid quantum-classical optimization algorithm which uses depicted in Fig. 4. This configuration is compatible with the
a unitary U(θ), characterized by the rotation angles θ, to NFT optimizer. For simulations, two-qubit gates are applied
prepare a parameterized quantum state |ψ(θ)⟩ = U(θ)|0⟩. in a “pairwise” fashion composed of two layers. In the first
The probability distribution of bitstrings, and the energies layer, qubit i is entangled with qubit (i + 1), for all even
correspondingtothem,sampledaftermeasurementofthestate values of i and in the second layer qubit i is entangled with
|0⟩ e−iθ0Y e−iθ3Y
|0⟩ e−iθ1Y Z e−iθ4Y
|0⟩ e−iθ2Y Z e−iθ5Y
repeat p times
Fig. 4: Two-local ansatz for a three-qubit case.
1.0
0.8
0.6
0.4
0.2
0.0
10 15 20 25 30 35 40
Number of qubits
ytilibaborp
sseccuS
60
50
40
30
20
10
0
10 15 20 25 30 35 40
Number of qubits
Fig. 6: Box plot of the optimality gap χ for problem sizes
between 10 and 40 qubits for CVaR-VQE runs on the MPS
simulator. The averaging scheme is the same as described in
Fig. 5.
the problems. The NFT optimizer evaluates n = 3 sets eval
of circuit parameters at each iteration. Each iteration consists
of n × N circuit calls. The optimization scheme is
eval shots
executed for N iterations until sufficient convergence of
iter
the CVaR values is observed, typically in the range of 100 to
200 iterations.
Fig.5:Boxplotdepictingthesuccesspercentagesforproblem
Because classical optimizers like NFT are prone to getting
sizes between 10 and 40 qubits for CVaR-VQE runs on
stuck in local minima, we run N = 10 independent
trial
the MPS simulator. For each problem size, simulations were
trials of CVaR-VQE for each problem instance. Each trial
conducted on ten distinct sequences. For each sequence, ten
run starts with a set of circuit parameters initialized randomly
independent CVaR-VQE trials were performed. The success
from a uniform distribution. We examine 10 different mRNA
percentageiscalculatedasthenumberofsuccessfulrunsover
sequences for each qubit size.
tentrials.Theaverageoverdifferentsequencesisindicatedby
Fig.5presentsaboxplotillustratingthesuccessprobability
blue triangles.
across various problem instances ranging from 10 and 40
qubits. A run is deemed successful if the bitstring with the
lowest energy, identical to the solution obtained by CPLEX,
qubit (i+1) for all odd values of i. For the hardware runs,
is found within the samples collected after measuring the
the two-qubit gates are applied between qubits belonging to
final quantum state |ψ(θ )⟩, where θ denotes the circuit
the actual hardware sub-graphs obtained after the hardware f f
parameters obtained after N iterations. The success prob-
characterization steps described in Section V-B1. iter
ability p is then defined as N /N , where N is
succ succ trial succ
V. RESULTS the number of successful trials. As evident from Fig. 5, the
CVaR-VQE algorithm demonstrates robust performance, with
A. Simulation runs
the average success probability in the range from 0.4 to 1.0.
We perform noise-free simulation of CVaR-VQE using the
In addition to the success probability, we also use the
matrix-product state (MPS) [57], [58] simulator provided in
optimalitygapasametrictoassessthequalityofthesolutions
Qiskit [46]. Since the computational cost of MPS-based sim-
obtainedfromCVaR-VQE.Theoptimalitygapχisdefinedas
ulator grows exponentially with entanglement [57], we limit
our simulations to shallow-depth circuits. We use two layers
|F(θ ) −F |
(p = 2) of the “two-local” ansatz as shown in Fig. 4. After χ= f low 0 ×100, (5)
|F |
examining several different values of the CVaR parameter α 0
(see Section IV-B) on relatively small problem sizes, we set it where F(θ ) represents the lowest value of the objective
f low
equal to 0.1 for all the simulation runs. The number of shots function, as defined in Eq. 3, found among the samples at the
N used in the simulations was 25 and 210 for the 10- end of N iterations. F denotes the value of the objective
shots iter 0
and 15-qubit problems, respectively, and 213 for the rest of function of the optimal solution found by CPLEX. χ is zero
TABLE I: Quantum device circuit and run attributes. Max
when the lowest energy found by the quantum optimization
iterations column indicates the maximum number of classi-
matches the one found by CPLEX, and greater than zero
cal iteration steps carried out by the NFT optimizer. The
otherwise. The smaller the deviation from zero, the better the
two-qubit Echoed Cross Resonant (ECR) gate is native to
quality of solutions obtained using CVaR-VQE.
ibm_brisbaneandibm_osaka.TheControl-Z(CZ)gate
Fig.6showsaboxplotoftheoptimalitygapχ ,averaged
avg is native to ibm_torino.
over 10 independent trials for each sequence, for qubits sizes
ranging from 10 to 40. The CVaR-VQE algorithm demon- ibm_brisbane/ibm_osaka ibm_torino
strates high performance, achieving an average optimality gap Qubits 26 40 50 80
of less than 20%. Circuitdepth 18 20 20 17
Maxiterations 104 320 600 –
Our results indicate that the performance of the algorithm ECR/CZCount 25 39 49 158
deteriorates as the problem size increases. For the constant-
depth circuits employed in this study, such a decline with
increasing problem size is anticipated. We hypothesize that
2) Error suppression and error mitigation: We use Dy-
the performance of CVaR-VQE can be further enhanced by
namical Decoupling (DD) [61], [62], [63], [64], which adds
increasing the number of layers in the ansatz, optimizing the
pulsesduringtheidletimesofthecircuitandhelpstosuppress
choice of α, or by incorporating other recent advancements
errors without any additional overhead. Additionally, we use
in CVaR [59]. We intend to investigate these aspects of the
a version of Matrix-free Measurement Mitigation (M3) [65]
algorithm in future research.
for readout error mitigation. Instead of working based off a
full assignment matrix of size 2n, M3 uses a much smaller
B. Hardware runs
subspaceformedbythenoisybitstrings,allowingittousually
handle large system efficiently [65]. In order to keep the
error mitigation overhead to “minimal”, we do not consider
|0⟩ H mitigation schemes like ZNE [23], [66], [67] and PEC [68]
in this work. Interestingly, as pointed out in Ref. [20], CVaR-
basedschemeswitharelativelysmallsamplingoverheadmay
|0⟩ H
still be able to track the noise-free expectation values [20].
We plan to explore these aspects of the CVaR algorithm in
SWAP
near future. Additionally, we set optimization level=3 under
Fig. 7: Quantum circuit for noise characterization.
theQiskit’sSamplerV1primitiveruntimeoptionstooptimize
circuits, which gives rise to transpiled circuits depths in the
In this Section, we provide details about the hardware runs range of 17-20 (see Table I).
carriedoutontheIBMEagleandHeronprocessors,including 3) HardwareResults: Onthehardware,westudy3mRNA
ibm_brisbane, ibm_osaka, and ibm_torino. sequences: (a) of length 25, represented by the problem
1) QubitSelection: Duetovariationsinthenoisecharacter- instance with 26 qubits, (b) of length 30, with 40 qubits, and
isticsacrossthephysicalchipofaquantumcomputer,mapping (c)oflength30,with50qubits.Werun8,8,and5independent
problem qubits to the hardware sub-graphs with smaller error CVaR-VQE trials for these instances, correspondingly, using
rates becomes crucial. To account for the cross-talk, similar α=0.25 and p=1.
in spirit to the layered fidelity benchmarking scheme [60], we The success probability of CVaR-VQE experiments on
examine fidelity across circuit layers containing disjoint set quantumhardwareisshowninFig.10.Itshouldbenotedthat
of qubit couplings. We implement a simple two-qubit circuit hardware experiments were conducted for only one sequence
(see Fig. 7) across couplings on each layer separately. Note, of each problem size; therefore, the success probability is
that the three layers with couplings of color blue, red, and computedasacumulativequantityover10trials.Asindicated
greeninFig.8coveralloftheneighboringinteractionspresent bytheplot,thesuccessprobabilityacrossdifferentqubitsizes
on the quantum device. Fidelity between ideal (00 bitstring) is above ∼ 0.38. The average optimality gap χ avg of CVaR-
and the distribution obtained from the quantum computer is VQE experiments on the hardware is shown in Fig. 11 and
indicated in Fig. 8 for the three layers. Next, in order to is below ∼34%. We set the number of shots for each circuit
select the set of “good” qubits, we consider a 1D-chain of to 28 for the 26-qubit problem and 213 for the larger problem
100 qubits, color coded as nodes in magenta in Fig. 8 [60]. instances.Themaximumnumberofiterationsrequiredtofind
Across the chain, we exclude qubits connected to couplings the optimal parameters on the hardware is reported in the
withfidelitiesnotmeetingasetthreshold.Forourrunsweset Table I. Some trials took fewer than the maximum number
this fidelity threshold to 0.85. If at the end of the process, we of iterations reported. The lowest energy sequences found are
still need more qubits, we extend our search to the rest of the shown in Fig. 9.
couplings not covered by the 100-qubit chain. We employ the For the 42-nucleotide, 80-qubit mRNA problem instance,
described process of qubit selection on all the devices once we do not carry out circuit parameter optimization on the
before starting the CVaR-VQE runs. hardware.Theoptimalparametersforthecircuitwereobtained
Coupling Index
(a)
32 36 51 70 74 89 108 112 126
13 31 50 69 88 107 125
12 17 30 49 55 68 87 93 106 124
11 29 48 67 86 105 123
10 28 35 47 66 73 85 104 111 122
9 27 46 65 84 103 121
8 16 26 45 54 64 83 92 102 120
7 25 44 63 82 101 119
6 24 34 43 62 72 81 100 110 118
5 23 42 61 80 99 117
4 15 22 41 53 60 79 91 98 116
3 21 40 59 78 97 115
2 20 33 39 58 71 77 96 109 114
1 19 38 57 76 95 113
0 14 18 37 52 56 75 90 94
(b)
Fig. 8: Panel a) shows the fidelity of different couplings across the three layers obtained during the hardware characterization
scheme detailed in Section V-B1. Panel b) shows the heavy-hex layout on ibm_brisbane and ibm_osaka with nodes on
the graph representing qubits and edges representing qubit couplings. The couplings belonging to the three layers are color
coded as blue, red, and, green, matching the fidelity plots in a). The purple colored nodes represent the 100 qubit chain
considered for selection of qubits on ibm_brisbane based on the fidelity values.
by running CVaR-VQE using the MPS simulator with the [76],[77],theheavy-hexlattice-basedcircuitscharacterizedby
α = 0.1 and p = 2. We use the optimized parameters to depths we employ in this study can be efficiently simulated
run circuits on the recently introduced ibm_torino device usingclassicalschemes.Ourresultspointtoanexcitingregime
featuring the 133-qubit tunable-coupler Heron processor [69]. where relatively large scale optimization-based problems, ap-
The resulting circuit has depth 17 and contains 158 two-qubit proaching utility-scale, can also be successfully executed on
CZ gates. We plot the distribution of measurement probabili- the noisy quantum hardware and can produce reliable results
ties of the sampled bitstrings with the corresponding energies using low-overhead error mitigation schemes.
in Fig. 12. The plot illustrates that the lowest energy bitstring
with energy value of -140, which matches the lowest energy C. Feasibility of hardware runs
bitstring found by CPLEX, is sampled with a probability
Unlike schemes targeting fault-tolerant quantum devices,
separated from the rest of the sampling probabilities by more
where quantum resource scaling arguments are relatively
than one order of magnitude. The measurement probability
straightforward to derive, establishing such arguments for
of the lowest bitstring is 0.045, hence, we need at most
variationalschemesismorechallenging[78].Ifplayersofan
1/0.045∼20samplestofindthelowestenergybitstring.The
ansatzareemployed,andeachcircuitlayerrequirestimeT to
optimal sequence found is shown in Fig. 12. Note, the offline
execute,thetotaltimetakenisapproximatelyN ×p×T,
training of circuit parameters using MPS simulator was only circuit
whereN representsthetotalnumberofcircuitsexecuted.
possiblebecauseoftheclassicalsimulatabilityofsuchcircuits. circuit
Giventhatthisstudyutilizesfixed-depthcircuitsandassuming
As pointed out in Refs. [21], [70], [71], [72], [73], [74], [75],
theexecutiontimeT remainsrelativelyinvariantwithproblem
1.0
0.8
0.6
0.4
0.2
0.0
26 40 50
Number of qubits
(a) 26 qubits, 25 nucleotides (b) 40 qubits, 30 nucleotides
(c) 50 qubits, 30 nucleotides
Fig. 9: Figures show actual optimal sequences based on the
lowest energy bitstring found during hardware runs.
size, the primary factor that necessitates closer examination is
N .
circuit
ThenumberofcircuitsN executedduringthescheme
circuit
isdeterminedbytheexpression1/p ×N ×N ×n ,
succ iter shots eval
where n represents the number of different sets of circuit eval
parameters evaluated at each iteration to suggest a new set
of circuit parameters θ1. To assess the feasibility of solving
structure prediction problems on a quantum computer, it is
crucial that none of the three factors {1/p ,N ,N },
succ iter shots
with n assumed constant, scale exponentially with the
eval
number of qubits. A concerning trend observed is the decline
in success probabilities p as the problem size increases.
succ
However, based on the results, it is difficult to definitively
concludethatthedeclineisexponential.Additionally,sincethe
1FortheNFToptimizer,n
eval
=3
ytilibaborp
sseccuS
Fig.10:Successprobabilityfordifferentnumberofqubitsfor
the hardware runs.
100
80
60
40
20
0
26 40 50
Number of qubits
Fig. 11: Box plot of the optimality gap χ for problem sizes in
the range of 10 to 40 qubits executed on quantum hardware.
Average χ over N independent runs is shown using blue
trial
triangles.
choiceofalgorithmdictatesthelowerboundonp ,thereis succ
potential for improvement through algorithmic advances. Un-
fortunately,beyondacertainregimeofclassicalsimulatability,
performancecomparisonofdifferentoptimizationschemescan
onlybecarriedoutinconjunctionwiththequantumhardware.
VI. CONCLUSIONSANDOUTLOOK
In this work, we examined the feasibility of solving sec-
ondary structure prediction-based optimization problems of
RNA sequences on a quantum computer. Our results demon-
strate that optimization schemes like CVaR-VQE running on
noisy quantum hardware with careful calibration and nominal
102
103
104
105
140 130 120 110 100 90 80 70
Energy
ytilibaborp
tnemerusaeM
more advanced optimization techniques.
Ourresultsindicatethatquantumcomputingoffersvaluable
insights into RNA folding mechanisms and enables rapid
identification of functional RNA structures. This work lays
thegroundworkforfutureresearchinquantum-assistedbioin-
formaticsandhighlightstheimportanceofdevelopingscalable
quantum algorithms for real-world biological problems.
ACKNOWLEDGMENTS
The authors would like to thank Yashrajsinh Jadeja and
Haining Lin for their valuable assistance in acquiring test
sequences and for their contributions to the prediction of
classical RNA secondary structures.
REFERENCES
[1] Q.VicensandJ.S.Kieft,“Thoughtsonhowtothink(andtalk)about
RNAstructure,”ProceedingsoftheNationalAcademyofSciences,vol.
119,no.17,p.e2112677119,2022.
[2] M. Metkar, C. S. Pepin, and M. J. Moore, “Tailor made: the art of
therapeutic mRNA design,” Nature Reviews Drug Discovery, vol. 23,
no.1,pp.67–83,2024.
[3] S. Qin, X. Tang, Y. Chen, K. Chen, N. Fan, W. Xiao, Q. Zheng,
G.Li,Y.Teng,M.Wuetal.,“mRNA-basedtherapeutics:powerfuland
versatile tools to combat diseases,” Signal Transduction and Targeted
Therapy,vol.7,no.1,p.166,2022.
[4] D. H. Mathews, W. N. Moss, and D. H. Turner, “Folding and finding
RNAsecondarystructure,”ColdSpringHarborPerspectivesinBiology,
vol.2,no.12,p.a003665,2010.
[5] R. B. Lyngsø and C. N. Pedersen, “RNA pseudoknot prediction in
energy-based models,” Journal of Computational Biology, vol. 7, no.
3-4,pp.409–427,2000.
[6] M. Zuker, “Computer prediction of RNA structure,” in Methods in
enzymology. Elsevier,1989,vol.180,pp.262–288.
[7] M.Zuker,“Mfoldwebserverfornucleicacidfoldingandhybridization
prediction,” Nucleic Acids Research, vol. 31, no. 13, pp. 3406–3415,
2003.
[8] P.Schuster,P.F.Stadler,andA.Renner,“RNAstructuresandfolding:
from conventional to new issues in structure predictions,” Current
OpinioninStructuralBiology,vol.7,no.2,pp.229–235,1997.
[9] M. Zuker, “On finding all suboptimal foldings of an RNA molecule,”
Science,vol.244,no.4900,pp.48–52,1989.
[10] E. P. Nawrocki and S. R. Eddy, “Infernal 1.1: 100-fold faster RNA
homology searches,” Bioinformatics, vol. 29, no. 22, pp. 2933–2935,
2013.
Fig. 12: (a) Measurement probability of sampled bitstrings [11] S.Zakov,Y.Goldberg,M.Elhadad,andM.Ziv-Ukelson,“Richparame-
terizationimprovesRNAstructureprediction,”JournalofComputational
plotted against bitstring energies for the 80-qubit problem.
Biology,vol.18,no.11,pp.1525–1542,2011.
Red bar indicates the probability corresponding to the lowest [12] M.Szikszai,M.Wise,A.Datta,M.Ward,andD.H.Mathews,“Deep
energy bitstring. (b) Optimal folded structure of the 42- learningmodelsforRNAsecondarystructureprediction(probably)do
not generalize across families,” Bioinformatics, vol. 38, no. 16, pp.
nucleotide, 80-qubit mRNA sequence based on the corre-
3892–3899,2022.
sponding lowest energy bitstring found by the hardware run. [13] H. Zhang, L. Zhang, D. H. Mathews, and L. Huang, “Linearpartition:
linear-timeapproximationofRNAfoldingpartitionfunctionandbase-
pairing probabilities,” Bioinformatics, vol. 36, no. Supplement 1, pp.
i258–i267,2020.
error-mitigation,canproducesufficientlyreliableoutputseven [14] D.M.Fox,C.M.MacDermaid,A.M.Schreij,M.Zwierzyna,andR.C.
for relatively large problem instances. However, the scala- Walker,“RNAfoldingusingquantumcomputers,”PLOSComputational
Biology,vol.18,no.4,p.e1010032,2022.
bility of these methods to larger problem sizes exceeding
[15] T. Zaborniak, J. Giraldo, H. Mu¨ller, H. Jabbari, and U. Stege, “A
100 qubits remains an open question that warrants further QUBO model of the RNA folding problem optimized by variational
investigation.Asquantumhardwarecontinuestoadvance,any hybridquantumannealing,”in2022IEEEInternationalConferenceon
Quantum Computing and Engineering (QCE). IEEE, 2022, pp. 174–
improvements in algorithmic techniques will further bolster
185.
prospects of using quantum computers for optimization. This [16] J.Jiang,Q.Yan,Y.Li,Y.Chai,M.Lu,Z.Cui,M.Dou,Q.Wang,Y.-C.
study, which primarily focused on the capabilities of utility- Wu,andG.-P.Guo,“PredictingRNAsecondarystructureonuniversal
quantumcomputer,”arXivpreprintarXiv:2305.09561,2023.
scalequantumhardware,onlyconsideredasubsetofavailable
[17] E. Farhi, J. Goldstone, and S. Gutmann, “A quantum approximate
quantum approaches; however, our future work will explore optimizationalgorithm,”arXivpreprintarXiv:1411.4028,2014.
[18] M. P. Harrigan, K. J. Sung, M. Neeley, K. J. Satzinger, F. Arute, [39] D. Draper, T. Gluick, and P. Schlax, “In RNA structure and function,
K.Arya,J.Atalaya,J.C.Bardin,R.Barends,S.Boixoetal.,“Quantum vol.298,”2000.
approximate optimization of non-planar graph problems on a planar [40] S. D. Gilbert, R. P. Rambo, D. Van Tyne, and R. T. Batey, “Structure
superconducting processor,” Nature Physics, vol. 17, no. 3, pp. 332– of the SAM-II riboswitch bound to s-adenosylmethionine,” Nature
336,2021. Structural&MolecularBiology,vol.15,no.2,pp.177–182,2008.
[19] Z. He, R. Shaydulin, S. Chakrabarti, D. Herman, C. Li, Y. Sun, and [41] D. J. Klein, T. E. Edwards, and A. R. Ferre´-D’Amare´, “Cocrystal
M.Pistoia,“AlignmentbetweeninitialstateandmixerimprovesQAOA structureofaclassIpreQ1riboswitchrevealsapseudoknotrecognizing
performance for constrained optimization,” npj Quantum Information, anessentialhypermodifiednucleobase,”NatureStructural&Molecular
vol.9,no.1,p.121,2023. Biology,vol.16,no.3,pp.343–344,2009.
[20] S. V. Barron, D. J. Egger, E. Pelofske, A. Ba¨rtschi, S. Eidenbenz, [42] R. C. Spitale, A. T. Torelli, J. Krucinska, V. Bandarian, and J. E.
M. Lehmkuehler, and S. Woerner, “Provable bounds for noise-free Wedekind,“ThestructuralbasisforrecognitionofthePreQ0metabolite
expectation values computed from noisy samples,” arXiv preprint byanunusuallysmallriboswitchaptamerdomain*,”JournalofBiolog-
arXiv:2312.00733,2023. icalChemistry,vol.284,no.17,pp.11012–11016,2009.
[21] S. H. Sack and D. J. Egger, “Large-scale quantum approximate opti- [43] D. H. Mathews, J. Sabina, M. Zuker, and D. H. Turner, “Expanded
mizationonnonplanargraphswithmachinelearningnoisemitigation,” sequencedependenceofthermodynamicparametersimprovesprediction
PhysicalReviewResearch,vol.6,no.1,p.013223,2024. ofRNAsecondarystructure,”JournalofMolecularBiology,vol.6,pp.
[22] A.Miessen,D.J.Egger,I.Tavernelli,andG.Mazzola,“Benchmarking 911–940,1999.
digitalquantumsimulationsandoptimizationabovehundredsofqubits [44] D.Gusfield,IntegerLinearProgramminginComputationalandSystems
using quantum critical dynamics,” arXiv preprint arXiv:2404.08053, Biology:AnEntry-LevelTextandCourse. CambridgeUniversityPress,
2024. 2019.
[23] Y.Kim,A.Eddins,S.Anand,K.X.Wei,E.VanDenBerg,S.Rosenblatt, [45] D. H. Turner and D. H. Mathews, “NNDB: the nearest neighbor
H. Nayfeh, Y. Wu, M. Zaletel, K. Temme et al., “Evidence for the parameter database for predicting stability of nucleic acid secondary
utilityofquantumcomputingbeforefaulttolerance,”Nature,vol.618, structure,”NucleicAcidsResearch,vol.38,no.suppl 1,p.D280–D282,
no.7965,pp.500–505,2023. 2009.
[24] P.K.Barkoutsos,G.Nannicini,A.Robert,I.Tavernelli,andS.Woerner, [46] Qiskit contributors, “Qiskit: An open-source framework for quantum
“ImprovingVariationalQuantumOptimizationusingCVaR,”Quantum, computing,”2023.
vol. 4, p. 256, Apr. 2020. doi: 10.22331/q-2020-04-20-256. [Online]. [47] A.Robert,P.K.Barkoutsos,S.Woerner,andI.Tavernelli,“Resource-
Available:https://doi.org/10.22331/q-2020-04-20-256 efficientquantumalgorithmforproteinfolding,”npjQuantumInforma-
[25] K.M.Nakanishi,K.Fujii,andS.Todo,“Sequentialminimaloptimiza- tion,vol.7,no.1,p.38,2021.
tion for quantum-classical hybrid algorithms,” Phys. Rev. Res., vol. 2, [48] I. Kolotouros and P. Wallden, “An evolving objective function
p.043158,Oct2020.doi:10.1103/PhysRevResearch.2.043158 for improved variational quantum optimisation,” arXiv preprint
[26] M.Oliv,A.Matic,T.Messerer,andJ.M.Lorenz,“Evaluatingtheimpact arXiv:2105.11766,2021.
of noise on the performance of the variational quantum eigensolver,” [49] W.Lavrijsen,A.Tudor,J.Mu¨ller,C.Iancu,andW.DeJong,“Classical
arXivpreprintarXiv:2209.12803,2022. optimizersfornoisyintermediate-scalequantumdevices,”in2020IEEE
[27] IBM,“IBMILOGCPLEXoptimizer International Conference on Quantum Computing and Engineering
https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex- (QCE). IEEE,2020,pp.267–277.
optimizer.” [50] A.Pellow-Jarman,S.McFarthing,I.Sinayskiy,A.Pillay,andF.Petruc-
[28] A.Mamuye,E.Merelli,andL.Tesei,“Agraphgrammarformodelling cione, “QAOA performance in noisy devices: the effect of classical
RNAfolding,”ElectronicProceedingsinTheoreticalComputerScience, optimizersandansatzdepth,”arXivpreprintarXiv:2307.10149,2023.
vol.231,pp.31–41,122016.doi:10.4204/EPTCS.231.3 [51] H. Singh, S. Majumder, and S. Mishra, “Benchmarking of different
[29] ViennaRNA,“ViennaRNAserverhttp://rna.tbi.univie.ac.at/.” optimizers in the variational quantum algorithms for applications in
[30] R. Durbin, S. R. Eddy, A. Krogh, and G. Mitchison, Biological Se- quantumchemistry,”TheJournalofChemicalPhysics,vol.159,no.4,
quence Analysis: Probabilistic Models of Proteins and Nucleic Acids. 2023.
Cambridgeuniversitypress,1998. [52] M.Ostaszewski,E.Grant,andM.Benedetti,“Structureoptimizationfor
[31] M. Andronescu, A. Condon, H. H. Hoos, D. H. Mathews, and K. P. parameterizedquantumcircuits,”Quantum,vol.5,p.391,2021.
Murphy, “Efficient parameter estimation for RNA secondary structure [53] W. Li, Y. Ge, S.-X. Zhang, Y.-Q. Chen, and S. Zhang, “Efficient and
prediction,”Bioinformatics,vol.23,no.13,pp.i19–i28,2007. robust parameter optimization of the unitary coupled-cluster ansatz,”
[32] J. Zuber, S. J. Schroeder, H. Sun, D. H. Turner, and D. H. Mathews, JournalofChemicalTheoryandComputation,2024.
“Nearest neighbor rules for RNA helix folding thermodynamics: im- [54] J. G. Vidal and D. O. Theis, “Calculus on parameterized quantum
provedendeffects,”NucleicAcidsResearch,vol.50,no.9,pp.5251– circuits,”arXivpreprintarXiv:1812.06323,2018.
5262,2022. [55] R.M.Parrish,J.T.Iosue,A.Ozaeta,andP.L.McMahon,“Ajacobidiag-
[33] M. W. Lewis, A. Verma, and T. T. Eckdahl, “Qfold: a new onalizationandandersonaccelerationalgorithmforvariationalquantum
modeling paradigm for the RNA folding problem,” Journal of algorithm parameter optimization,” arXiv preprint arXiv:1904.03206,
Heuristics, vol. 27, pp. 695–717, 2021. [Online]. Available: https: 2019.
//doi.org/10.1007/s10732-021-09471-3 [56] D.Wierichs,J.Izaac,C.Wang,andC.Y.-Y.Lin,“Generalparameter-
[34] J.-H. Chen, S.-Y. Le, and J. V. Maizel, “Prediction of common shiftrulesforquantumgradients,”Quantum,vol.6,p.677,2022.
secondary structures of RNAs: a genetic algorithm approach,” Nucleic [57] G. Vidal, “Efficient classical simulation of slightly entangled quantum
Acids Research, vol. 28, no. 4, pp. 991–999, 02 2000. doi: computations,”PhysicalReviewLetters,vol.91,no.14,p.147902,2003.
10.1093/nar/28.4.991. [Online]. Available: https://doi.org/10.1093/nar/ [58] U.Schollwo¨ck,“Thedensity-matrixrenormalizationgroupintheageof
28.4.991 matrixproductstates,”AnnalsofPhysics,vol.326,no.1,pp.96–192,
[35] Q. Zhao, Z. Zhao, X. Fan, Z. Yuan, Q. Mao, and Y. Yao, “Review 2011.
ofmachinelearningmethodsforRNAsecondarystructureprediction,” [59] I.KolotourosandP.Wallden,“Evolvingobjectivefunctionforimproved
PLoSComputBiol,vol.17,no.8:e1009291,2021.[Online].Available: variational quantum optimization,” Physical Review Research, vol. 4,
https://doi.org/10.1371/journal.pcbi.1009291 no.2,p.023225,2022.
[36] D. H. Seetin, Matthew G.and Mathews, RNA Structure Prediction: [60] D. C. McKay, I. Hincks, E. J. Pritchett, M. Carroll, L. C. Govia, and
An Overview of Methods. Totowa, NJ: Humana Press, 2012, S.T.Merkel,“Benchmarkingquantumprocessorperformanceatscale,”
pp. 99–122. ISBN 978-1-61779-949-5. [Online]. Available: https: arXivpreprintarXiv:2311.05933,2023.
//doi.org/10.1007/978-1-61779-949-5 8 [61] L.ViolaandS.Lloyd,“Dynamicalsuppressionofdecoherenceintwo-
[37] S. R. Eddy, “How do RNA folding algorithms work?” Nature state quantum systems,” Physical Review A, vol. 58, no. 4, p. 2733,
Biotechnology, pp. 1457–1458, 2004. doi: 10.1038/nbt1104-1457. 1998.
[Online].Available:https://doi.org/10.1038/nbt1104-1457 [62] P.Zanardi,“Symmetrizingevolutions,”PhysicsLettersA,vol.258,no.
[38] M. Zuker and P. Stiegler, “Optimal computer folding of large RNA 2-3,pp.77–82,1999.
sequences using thermodynamics and auxiliary information,” Nucleic [63] D.VitaliandP.Tombesi,“Usingparitykicksfordecoherencecontrol,”
AcidsResearch,vol.9,no.1,pp.133–148,1981. PhysicalReviewA,vol.59,no.6,p.4178,1999.
[64] N. Ezzell, B. Pokharel, L. Tewala, G. Quiroz, and D. A. Lidar, “Dy-
namicaldecouplingforsuperconductingqubits:aperformancesurvey,”
PhysicalReviewApplied,vol.20,no.6,p.064027,2023.
[65] P. D. Nation, H. Kang, N. Sundaresan, and J. M. Gambetta, “Scalable
mitigationofmeasurementerrorsonquantumcomputers,”PRXQuan-
tum,vol.2,no.4,p.040326,2021.
[66] A. Kandala, K. Temme, A. D. Co´rcoles, A. Mezzacapo, J. M. Chow,
andJ.M.Gambetta,“Errormitigationextendsthecomputationalreach
ofanoisyquantumprocessor,”Nature,vol.567,no.7749,pp.491–495,
2019.
[67] Y. Li and S. C. Benjamin, “Efficient variational quantum simulator
incorporating active error minimization,” Physical Review X, vol. 7,
no.2,p.021050,2017.
[68] E. Van Den Berg, Z. K. Minev, A. Kandala, and K. Temme, “Proba-
bilistic error cancellation with sparse Pauli–Lindblad models on noisy
quantum processors,” Nature Physics, vol. 19, no. 8, pp. 1116–1121,
2023.
[69] IBM Newsroom. IBM debuts next-generation quantum processor
& IBM quantum system two, extends roadmap to advance era
of quantum utility. [Online]. Available: https://newsroom.ibm.com/
2023-12-04-IBM-Debuts-Next-Generation-Quantum-Processor-IBM-Quantum-System-Two,
-Extends-Roadmap-to-Advance-Era-of-Quantum-Utility
[70] S.Patra,S.S.Jahromi,S.Singh,andR.Oru´s,“Efficienttensornetwork
simulation of IBM’s largest quantum processors,” Physical Review
Research,vol.6,no.1,p.013326,2024.
[71] J. Tindall, M. Fishman, E. M. Stoudenmire, and D. Sels, “Efficient
tensor network simulation of IBM’s Eagle kicked Ising experiment,”
PRXQuantum,vol.5,no.1,p.010308,2024.
[72] H.-J.Liao,K.Wang,Z.-S.Zhou,P.Zhang,andT.Xiang,“Simulationof
IBM’skickedIsingexperimentwithprojectedentangledpairoperator,”
arXivpreprintarXiv:2308.03082,2023.
[73] S. Anand, K. Temme, A. Kandala, and M. Zaletel, “Classical bench-
marking of zero noise extrapolation beyond the exactly-verifiable
regime,”arXivpreprintarXiv:2306.17839,2023.
[74] T. Begusˇic´ and G. K. Chan, “Fast classical simulation of evidence for
theutilityofquantumcomputingbeforefaulttolerance,”arXivpreprint
arXiv:2306.16372,2023.
[75] T. Begusˇic´, J. Gray, and G. K.-L. Chan, “Fast and converged classical
simulations of evidence for the utility of quantum computing before
faulttolerance,”ScienceAdvances,vol.10,no.3,p.eadk4321,2024.
[76] M. S. Rudolph, E. Fontana, Z. Holmes, and L. Cincio, “Classical
surrogate simulation of quantum systems with lowesa,” arXiv preprint
arXiv:2308.09109,2023.
[77] Y.Shao,F.Wei,S.Cheng,andZ.Liu,“Simulatingquantummeanvalues
innoisyvariationalquantumalgorithms:Apolynomial-scaleapproach,”
arXivpreprintarXiv:2306.05804,2023.
[78] G. Scriva, N. Astrakhantsev, S. Pilati, and G. Mazzola, “Challenges
of variational quantum optimization with measurement shot noise,”
PhysicalReviewA,vol.109,no.3,p.032408,2024.
