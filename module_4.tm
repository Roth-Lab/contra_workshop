<TeXmacs|1.99.8>

<style|<tuple|tmbook|british>>

<\body>
  <section|Phylogenetic analysis>

  <subsection|Overview>

  In this module we will consider the problem of reconstructing the
  evolutionary tree or <with|font-shape|italic|phylogeny> that relates cancer
  clones. First we will review the types of data we have available to address
  this problem. Next we do a brief review of probabilistic models in
  phylogenetics. Finally, we develop a probabilistic model to infer
  phylogenies from bulk. This model has a few unique twists that are tailored
  to cancer genomics data.

  <subsection|Data types>

  The model we develop in this module will use bulk sequence data. Thus we
  focus the bulk of our discussion on this. However, single cell sequencing
  is becoming increasingly common and we offer some opinions on the potential
  challenges here.

  <subsubsection|Bulk sequencing>

  The problem of building clonal phylogenies is a very active area of
  research in computational cancer biology. To date most work has focused on
  using bulk sequencing to address this issue. The major challenge of using
  bulk sequencing is that we need to deal with fact samples are mixtures of
  cell populations. There are two different approaches that are widely
  employed to do this.

  The first approach is to perform deconvolution using methods such as the
  one we developed in module 2. The key assumption that all methods of this
  type make is that more prevalent mutations should originate higher up the
  tree. Within this category there are two different strategies. The first
  strategy is to perform the deconvolution and tree building simultaneously.
  The second strategy is to first perform the deconvolution and then
  reconstruct trees based on cellular prevalence. In general the first
  strategy is more computationally demanding. The key benefit is that the
  tree structure can inform the deconvolution step.

  The second general approach is to use classical phylogenetic methods. In
  standard phylogenetic problems we observe a set of species with an
  associated matrix of observations. The values in this matrix can either be
  binary, for example the presence or absence of a morphological feature.
  Alternatively they can be nucleotide or amino acids if we have aligned
  sequence data. In the context of cancer it is common to use the
  presence/absence of a mutation within a sample as the feature matrix. Here
  species correspond to samples. The obvious problem with this approach is
  that bulk samples represent mixtures of clonal populations. Thus there is
  no guarantee that the presence of a pair of mutations in sample implies
  they exist within the same cell. Because of this one need to be careful
  when interpreting <with|font-shape|italic|sample trees>. Naively
  constructed trees based on the approach outlined above simply represent
  similarities between samples, not an evolutionary relationship. Later in
  this module we will see one way to address this issue.

  One might wonder why anyone uses the second approach? Probably the dominant
  reason is that it is relatively simple, and does communicate some concept
  of similarity between samples. A more fundamental, but perhaps under
  appreciated reason is that we can account for mutation loss in this setup.
  In contrast, the deconvolution approach hinges on the fact that higher
  cellular prevalence mutations appear higher in the tree. This assumption is
  invalidated if mutations are lost. In practice many methods are robust to
  this issue. The basic strategy is to treat some mutations as outliers and
  not fit them to the tree. While practical, it is somewhat inelegant as we
  are using a model we believe to be wrong.

  <subsubsection|Single cell sequencing>

  The emergence of single cell sequencing provides hope we can avoid the
  clonal mixture issue that afflicts bulk sequencing. In principle single
  cell data is ideal for classical phylogenetic methods since cells are
  equivalent to species or individuals within a population. The major
  problem, as with all things single cell, is noise. To build a phylogeny we
  need a set of features for the observation matrix. The question is what
  could these be?\ 

  We could use SNVs, but low coverage whole genome sequencing is unlikely to
  have even a single read covering an SNV in a cell. This leads to a very
  sparse observation matrix. We could potentially use CNVs which are easily
  detectable from binned read count data even in low coverage sequencing.
  However, modelling CNV evolution is a formidable problem since we can no
  longer make the infinite sites assumption. The main issue is that we cannot
  treat the copy number of a bin as a feature, as independent events may have
  lead to the same observed copy number. One promising approach is to use the
  change-points associated with CNVs. These are less likely to violate the
  infinite sites assumption.\ 

  Scalability is another issue for single cell data. Typical phylogenetic
  methods become computationally demanding with 100s of species. We will soon
  have data sets with thousands of cells if not more. One option is to
  abandon probabilistic models, and use faster methods based on distance
  metrics. This is not without issue as the distance matrix may be noisy. An
  alternative is to use more advanced inference techniques. Sequential Monte
  Carlo (SMC) which we will discuss in the next module is one promising
  approach for Bayesian phylogenetics. There is also some promising work
  using variational inference to improve MCMC samples in a post-processing
  step.

  <subsubsection|Summary>

  Phylogenetic reconstruction from bulk sequencing data is challenging,
  primarily due to the need to deal with sample mixtures. Phylogenetics using
  single cell sequencing in contrast is conceptually simple, but challenging
  due to noise and scale. In addition, there is a great deal more bulk data
  currently available than single cell. Thus bulk phylogenetic methods are
  likely to remain important for the foreseeable future. With that in mind we
  will turn to the problem of reconstructing phylogenies from bulk sequence
  data.

  <subsection|Probabilistic phylogenetic models>

  Here we provide a brief review of probabilistic models for phylogenetics.
  First, we introduce the basic problem setup. Next we look at how define a
  probabilistic model and compute the probability of the data efficiently. We
  finish with a brief discussion of Bayesian models and the computational
  challenges.

  <subsubsection|Tree topology>

  The basic problem of phylogenetics is as follows. Given an observation
  matrix <math|X> with <math|M> rows corresponding to species and <math|N>
  columns corresponding to features (characters), we would like to infer the
  evolutionary tree relating the species. For simplicity we will assume this
  is a binary matrix i.e. <math|X\<in\><around*|{|0,1|}><rsup|M\<times\>N>>.
  Let the topology of the tree be denoted by <math|\<tau\>=<around*|(|E,V|)>>
  where <math|E> is the set of edges and <math|V> the set of vertices. In
  many problems there will also be branch lengths, <math|\<Lambda\>>,
  associated with each edge, though not in all cases. If we know the
  ancestral sequence or can identify a distantly related species to use as an
  outgroup we can root a tree. While this is not easy to do when looking at
  species evolution, in cancer it typically straightforward as we know the
  genotype of the normal progenitor cell (root) of the cancer. Thus we will
  focus strictly on rooted trees in this module. In a rooted tree we have a
  well defined notion of parent child relationships, and thus which nodes are
  ancestral.

  <subsubsection|Transition probabilities>

  The first step to defining a probabilistic phylogenetic model is to
  determine the probability of transitioning from one state to the other, in
  the simple case from a 0 to 1 or vice-verse. If we have branch lengths the
  probability of this transition will need to depend on the length of the
  branch. The usual approach is to define a rate matrix <math|Q> where
  <math|Q<rsub|i j>> is the instantaneous rate of transition from state
  <math|i> to <math|j>. We can then define a transition matrix as
  <math|P=exp<around*|(|Q t|)>> where <math|t> is time or branch length. Here
  <math|P<rsub|i j>> is the probability of a transition from state <math|i>
  to <math|j> along the branch with length <math|t>. If we do not have branch
  lengths we can define <math|P> directly.

  <\remark>
    Identifiability crops up again here. If we rescale the matrix <math|Q> by
    a constant and divide the branche length by the same value <math|P>
    remains unchaged. Biologically rescaling <math|Q> corresponds to altering
    the mutation rate. Informally this means we would expect to see the same
    number of mutations if the rate is high and time interval short as we
    would if the rate was low and time interval long.
  </remark>

  <subsubsection|Tree probability>

  Given a tree and transition probabilities the question is then how to
  define the probability of the data? For the moment assume we observe only a
  single character so that our data is a vector <math|\<b-x\>>. We would like
  to compute <math|p<around*|(|\<b-x\>\|\<tau\>,P|)>>, that is the
  probability of the data given the tree <math|\<tau\>> and the transition
  matrix <math|P>. We will assume that we know not only the character
  sequence of the leafs (species), but also of the unobserved internal nodes
  of the tree. Let

  <\itemize>
    <item><math|x<rsub|v>> denotye the value of the characet at leaf node
    <math|v>

    <item><math|y<rsub|v>> denote the value of the character at internal node
    <math|v>

    <item><math|L<around*|(|\<tau\>|)>> denote the leaf nodes of
    <math|\<tau\>>

    <item><math|I<around*|(|\<tau\>|)>=V\\L<around*|(|v|)><big|cup><around*|{|r|}>>
    denote the set of internal nodes that are not the root

    <item><math|\<tau\><around*|(|v|)>> denote the subtree rooted at node
    <math|v>

    <item><math|\<gamma\><around*|(|v|)>> denote the set of children of
    <math|v>

    <item><math|r> denote the root node of the tree
  </itemize>

  Then we have\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>\|\<tau\>,P,\<b-y\>|)>>|<cell|=>|<cell|<big|prod><rsub|v\<in\>L<around*|(|v|)>>P<rsub|y<rsub|\<rho\><around*|(|v|)>>
    x<rsub|v>> <big|prod><rsub|v\<in\>I<around*|(|\<tau\>|)>>P<rsub|y<rsub|\<rho\><around*|(|v|)>>
    y<rsub|v>>>>|<row|<cell|>|<cell|=>|<cell|<big|prod><rsub|v\<in\>\<gamma\><around*|(|r|)>>P<rsub|y<rsub|r>
    \ y<rsub|v>> p<around*|(|\<b-x\>\|\<tau\><around*|(|v|)>,P,\<b-y\>|)>>>>>
  </eqnarray>

  where we have developed a recursive definition in the second line. The
  recursion is vital to efficient computation as we will see later. We of
  course do not know the values of the unobserved hidden nodes so we would
  like to marginalise them, that is sum over all possible states of
  <math|\<b-y\>> to obtain

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>\|\<tau\>,P|)>>|<cell|=>|<cell|<big|sum><rsub|\<b-y\>>p<around*|(|\<b-x\>\|\<tau\>,P,\<b-y\>|)>>>>>
  </eqnarray>

  This can be done efficiently using the Felsenstein pruning algorithm, which
  is another example of dynamic programming. The key quantity is\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|\<alpha\><rsub|v><around*|(|s|)>>|<cell|=>|<cell|<big|prod><rsub|u\<in\>\<gamma\><around*|(|v|)>><big|sum><rsub|y<rsub|u>>P<rsub|s
    y<rsub|u>> \<alpha\><rsub|u><around*|(|y<rsub|u>|)>>>>>
  </eqnarray>

  and we define\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|\<alpha\><rsub|v><around*|(|s|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|1>|<cell|<text|if>>|<cell|s=x<rsub|v>>>|<row|<cell|0>|<cell|<text|if>>|<cell|s\<neq\>x<rsub|v>>>>>>>>>>
  </eqnarray>

  for leaf nodes. In words, the pruning algorithm recursively computes
  <math|\<alpha\><rsub|s><around*|(|v|)>> from the leaf nodes upwards. In
  each round we sum over all possible states for the current node and
  <math|\<alpha\><around*|(|s|)>> then take the product of all the child
  terms.

  <\remark>
    We can also modify the above recursion to compute the most probable set
    of internal states for the tree in the same way as we use the Viterbi
    algorithm for HMMs.
  </remark>

  <\remark>
    As we will see later it is not necessary that we observe the character
    values with complete certainty. The 0/1 condition for the leaf nodes can
    be replaced with probabilities.
  </remark>

  <subsubsection|Bayesian phylogenetics>

  Now that we can define the probability for the data, we need to specify
  priors to complete the Bayesian model. This includes a prior for the tree
  topology <math|p<around*|(|\<tau\>|)>>, the branch lengths if they are
  included <math|p<around*|(|\<Lambda\>|)>> and any additional parameters
  which govern the transition probabilities <math|p<around*|(|\<theta\>|)>>.
  The joint probability is then given by

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|X,\<tau\>,\<Lambda\>,\<theta\>|)>>|<cell|=>|<cell|p<around*|(|X\|\<tau\>,\<Lambda\>,\<theta\>|)>p<around*|(|\<tau\>|)>p<around*|(|\<Lambda\>|)>p<around*|(|\<theta\>|)>>>>>
  </eqnarray>

  \;

  The question now is how to fit the model? The most challenging aspect is
  inferring the tree topology <math|\<tau\>>. Even computing the MAP estimate
  is hard in most cases as it requires searching over the space of all
  possible trees. For rooted binary trees with <math|n> leafs there are
  <math|<frac|<around*|(|2n-3|)>!|2<rsup|n-2>
  <around*|(|n-2|)>!>\<approx\>n!> possible trees. This renders exhaustive
  search impossible for all but small <math|n>. This is one case where full
  Bayesian inference via MCMC is not much slower than other approaches. The
  typical way to perform inference is to use MH moves for <math|\<tau\>>.
  Proposing a new random tree rarely works well, so local moves such pruning
  a subtree and regrafting to another part of the tree are used in practice.
  The other parameters are typically easier to sample, and an array of MCMC
  \ methods such as MH or Hamiltoniain Monte Carlo can be applied.

  <\remark>
    Implementing MCMC for phylogenetics is a bit painful due to the book
    keeping required for moves in the tree space. Luckily many good tools
    such as Mr. Bayes and BEAST exist. These software packages provide a
    great deal of flexibility in defining models, and are worth looking into
    before developing your own tools.
  </remark>

  <subsection|Probabilistic model for mutation loss>

  We now turn to developing our own phylogenetic method in this section. We
  will address the problem of building trees from bulk sequencing. We first
  define the problem, then discuss a suitable model and finally look at some
  results from real world data.

  <subsubsection|Problem statement>

  We consider the problem of building phylogenetic trees from bulk sequence
  data. We will assume we have multiple tissue samples collected from the
  same patient. We will also assume CNV analysis has been performed and
  estimates of tumour content are available. Our input data will allelic
  counts from SNVs. We will use this data to compute an observation matrix
  where the entries are the probability a mutation is
  <with|font-shape|italic|clonally present>.

  We will assume that there is ongoing genomic instability in the cancer that
  deletes mutations. This will require us to define a non-standard model
  which is closely related to the Stochastic Dollo process first used in
  linguistics. This model assumes that mutations originate only once on the
  tree, but can be subsequently lost any number of times.

  <subsubsection|Model description>

  We would like to infer a rooted bifuricating tree relating a set of
  <math|M> samples where we observe <math|N> mutations. We assume that
  mutations originate only once on the tree at node <math|w>, which we cannot
  observe. Once a mutation has originated it is propogated to its children
  with probability <math|<around*|(|1-\<pi\><rsub|l>|)>> or lost with
  probability <math|\<pi\><rsub|l>>. Once the mutation is lost it remains
  lost in all descendants. There are no branch lengths in the model so the
  only parameter governing the transition probabilities is
  <math|\<pi\><rsub|l>>. The main assumptions are

  <\enumerate-numeric>
    <item>Mutations originate at most once on the tree

    <item>Mutations can be lost after they are acquired

    <item>Mutations evolve indepdently i.e. our tree probability decomposes
    as the product of mutations
  </enumerate-numeric>

  \;

  One unusual feature of the model is that we will not consider the
  observation data to be perfect. Rather we will assume there is a
  probability that a sample contains the mutations and denote this
  <math|p<around*|(|z<rsub|i>\|\<cdot\>|)>>. Our data matrix is then a set of
  probabilities that a mutation is present in a sample. We can revert to the
  simple case by setting <math|p<around*|(|z<rsub|i>\|\<cdot\>|)>> to one for
  the observed value and 0 for all others. For example we can threshold on
  the number of reads supporting the mutation.

  There are two reasons for using probabilities in the observation data.
  First, we assume the data from whole genome sequencing so read depth is on
  the order of 30x-100x. Depending on tumour content this means there is a
  reasonable chance we will fail to detect the mutation when it is present.
  We will see when we compute the likelihood a mutation is present we
  explicitly correct for tumour content. The second reason is to address the
  problem of samples representing mixtures of cells. As discussed earlier
  this raises issues when intepreting sample trees. A solution to this
  problem would be to only build trees with mutations that are clonal in the
  sample. This would guarantee that the mutations co-occur within the same
  cell. We do not have this information in practice, but using the same
  approach we developed in module 2 we can compute a probability for this.

  <subsubsection|Probability of clonal presence>

  To generate our data matrix, <math|X>, for phylogenetic reconstruction we
  need to compute the pobability a mutation is clonally present. We define
  this to be the probability a mutation has cellular prevalence 1.0 in the
  sample. To compute this probability we adopt an approach similar to module
  <math|2>. Let <math|c<rsub|b>> denote the number of mutated copies of a
  locus and <math|c<rsub|t>> the total number of copies. Let <math|t> be the
  tumour content of the sample. Then the probability of observing a read with
  the mutation, if it had cellular prevalence 1.0, is\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|r>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|<frac|c<rsub|b>
    t|2 <around*|(|1-t|)>+c<rsub|t> t>>|<cell|<text|if>>|<cell|c<rsub|b>\<gtr\>0>>|<row|<cell|\<varepsilon\>>|<cell|<text|if>>|<cell|c<rsub|b>=0>>>>>>>>>
  </eqnarray>

  so the probability of observing <math|b> reads with variant out of <math|d>
  total is\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|b\|d,c<rsub|b>,c<rsub|t>,t|)>>|<cell|=>|<cell|<text|Binomial><around*|(|b\|d,r|)>>>>>
  </eqnarray>

  Since we do not know the number of copies of the variant, we sum over all
  possible values of <math|c<rsub|b>\<in\><around*|{|1,\<ldots\>,c<rsub|t>|}>>
  to compute the probability of presence. We assume a uniform prior for
  <math|c<rsub|b>> so all values have equal prior weight. For the probability
  of absence we set <math|c<rsub|b>=0>.

  <subsubsection|Tree notation>

  <\itemize-dot>
    <item><math|V<around*|(|\<tau\>|)>> denote the vertices of the tree
    <math|\<tau\>>

    <item><math|L<around*|(|\<tau\>|)>> denote the leaves of the tree

    <item><math|D<around*|(|i|)>> denote the nodes descendant from node
    <math|i>

    <item><math|L<around*|(|i|)>> denote the leaves descendant from node
    <math|i>

    <item><math|C<around*|(|i|)>> denote the childrend of node <math|i>

    <item><math|\<rho\><around*|(|i|)>> denote the parento fo node <math|i>

    <item><math|A<around*|(|i|)>> denote the ancestors of node <math|i> i.e.
    all nodes on the path from <math|i> to the root

    <item><math|w> denote the tree node at which a mutation originated

    <item><math|z<rsub|i>> be an indicator if node <math|i> has the mutation

    <item><math|\<pi\><rsub|l>> be the probability of losing a mutation along
    an edge

    <item><math|p<around*|(|z<rsub|i>\|\<cdot\>|)>> be the likelihood a
    variant is present
  </itemize-dot>

  <subsubsection|Tree probability>

  We will develop a slightly non-standard recursion for computing tree
  probabilites. The reason for this is that the single origin constraint adds
  dependencies between branches in the tree. If a mutation originated in one
  branch it cannot originate again in another branch. Conditioned on the
  originating branch, <math|(\<rho\>(w), w)>, the losses are Markovian on the
  sub-tree <math|D(w)>. Informally this means once we know where a mutation
  orginiates evolution is independent among the nodes. Note that once a
  mutation is lost it will be lost in all descendants as well due to the
  single origin assumption.

  Mutations can originate at any point along the branch of the tree. To
  simplify the discussion we say a mutation originates at a node if it occurs
  at any point along the branch between the node at its parent. We use the
  same convention for lost mutations. As we will see later we can then infer
  at which node (clone) a mutation originated, and which clones subsequently
  lost the mutation.

  To compute the probability of the tree we imagine picking a node <math|w>
  as the node at which a mutation originated. Then we would like to compute
  <math|p<around*|(|\<b-x\>\|\<tau\>,w|)>> the probability of the data given
  the tree <math|\<tau\>> and origin node <math|w>. To do this we introduce
  the function <math|Q<around*|(|i,\<tau\>|)>>, the probability the mutation
  is present at node <math|i> given all possible combinations of losses on
  the sub-tree rooted at <math|i>. This can be compute recursively as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|Q<around*|(|j,\<tau\>|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|\<pi\><rsub|l>
    p<around*|(|z<rsub|j>=0<mid|\|>\<cdot\>|)>+<around*|(|1-\<pi\><rsub|l>|)>
    p<around*|(|z<rsub|j>=1<mid|\|>\<cdot\>|)>>|<cell|<text|if>>|<cell|j\<in\>L<around*|(|T|)>>>|<row|<cell|\<pi\><rsub|l>
    <big|prod><rsub|i\<in\>L<around*|(|j|)>>p<around*|(|z<rsub|i>=0\|\<cdot\>|)>+<around*|(|1-\<pi\><rsub|l>|)>
    <big|prod><rsub|i\<in\>C<around*|(|j|)>>Q<around*|(|i,\<tau\>|)>>|<cell|<text|if>>|<cell|j\<nin\>L<around*|(|T|)>>>>>>>>>>
  </eqnarray>

  The first line is the initial condition for the leaf nodes. The term
  <math|\<pi\><rsub|l> p<around*|(|z<rsub|j>=0<mid|\|>\<cdot\>|)>> accounts
  for losses on the branch above the node. The second term
  <math|<around*|(|1-\<pi\><rsub|l>|)> p<around*|(|z<rsub|j>=1<mid|\|>\<cdot\>|)>>
  accounts for the probability of not being lost and observed. The second
  line is the internal node recursion. The first term, <math|\<pi\><rsub|l>
  <big|prod><rsub|i\<in\>L<around*|(|j|)>>p<around*|(|z<rsub|i>=0\|\<cdot\>|)>>,
  again accounts for a loss but also enforces the loss for all children. The
  second term, <math|<around*|(|1-\<pi\><rsub|l>|)>
  <big|prod><rsub|i\<in\>C<around*|(|j|)>>Q<around*|(|i,\<tau\>|)>>, accounts
  for the event that the mutation is not lost and then considers all possible
  loss patterns on the child subtrees.\ 

  With <math|Q<around*|(|j,\<tau\>|)>> defined we then have

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>\|\<tau\>,w|)>>|<cell|=>|<cell|Q<around*|(|j=w,\<tau\>|)>
    <big|prod><rsub|i\<in\>L<around*|(|\<tau\>|)>\\L<around*|(|w|)>>p<around*|(|z<rsub|i>=0\|\<cdot\>|)>>>>>
  </eqnarray>

  which decomposes into a term for the subtree the mutation originates on,
  and term for all other nodes. Ultimately we want to compute
  <math|p<around*|(|\<b-x\>\|\<tau\>|)>> which can be obtained by
  marginalizing <math|w> over all nodes in the tree. In the absence of
  additionall information we assume a uniform prior for <math|w> i.e.
  <math|p<around*|(|w|)>=<frac|1|<around*|\||V<around*|(|T|)>|\|>>>. Thus we
  have

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>\|\<tau\>|)>>|<cell|=>|<cell|<big|sum><rsub|w\<in\>V<around*|(|T|)>>p<around*|(|\<b-x\>\|\<tau\>,w|)>
    p<around*|(|w|)>>>>>
  </eqnarray>

  To obtain the probability for all the data we use the assumption mutations
  are independent. Thus the probability of the data is given by

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|X\|\<tau\>|)>>|<cell|=>|<cell|<big|prod><rsub|n=1><rsup|N>p<around*|(|\<b-x\><rsub|n>\|\<tau\>|)>>>>>
  </eqnarray>

  where <math|\<b-x\><rsub|n>> is the data for the <math|n<rsup|th>>
  mutation.

  <subsubsection|Inference>

  We need to infer two parameters in this model: the tree topology
  <math|\<tau\>> and the probability of mutation loss <math|\<pi\><rsub|l>>.
  Fitting this model is potentially challenging because of the need to
  explore the space of trees. As stated earlier MCMC methods are usually
  favoured for addressing this problem. However, it is rare to perform
  multi-region sequencing on more than ten samples. This means that the
  number of trees is at most on the order of millions. In addition the
  computation of the tree likelihood can be performed in parallel for each
  tree topology. So we have a rare example in phylogenetics where performing
  MAP estimation is reasonably easy. To do this we simply enumerate all trees
  and optimize the mutation loss parameter for each tree. We then take the
  one with the highest joint probability as our solution.

  <\remark>
    Since we first developed this model we have generated larger datasets.
    Even with 12 samples the MAP estimation approach is no longer
    computationally feasible. An area of future work is to implement a proper
    MCMC sampler for this method.
  </remark>

  <subsubsection|Inferring origin, presence and loss of mutations>

  Given the MAP estimates for the tree <math|<wide|\<tau\>|^>> and the
  probability of loss <math|<wide|\<pi\>|^><rsub|l>> we can compute the most
  likely node where a mutation originated. We can also identify which at
  which nodes the mutation is present and at which the mutation was lost. The
  strategy is to maximize the following quantity

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|w,z\|x,<wide|\<tau\>|^>,<wide|\<pi\>|^><rsub|l>|)>>|<cell|=>|<cell|p<around*|(|x\|w,z,<wide|\<tau\>|^>,<wide|\<pi\>|^><rsub|l>|)>p<around*|(|z\|w|)>p<around*|(|w|)>>>>>
  </eqnarray>

  To achieve this we modify the previous recusion for
  <math|Q<around*|(|i,\<tau\>|)>> a follows

  <\eqnarray>
    <tformat|<table|<row|<cell|<wide|Q|^><around*|(|j,\<tau\>|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|max<around*|{|\<pi\><rsub|l>
    p<around*|(|z<rsub|j>=0\|\<cdot\>|)>,<around*|(|1-\<pi\><rsub|l>|)>
    p<around*|(|z<rsub|j>=1\|\<cdot\>|)>|}>>>|<row|<cell|max<around*|{|\<pi\><rsub|l>
    <big|prod><rsub|i\<in\>L<around*|(|j|)>>p<around*|(|z<rsub|i>=0\|\<cdot\>|)>,<around*|(|1-\<pi\><rsub|l>|)>
    <big|prod><rsub|i\<in\>C<around*|(|j|)>><wide|Q|^><around*|(|i,\<tau\>|)>|}>>>>>>>>>>
  </eqnarray>

  In other words we replace summation with maximisation in the recursion.
  This is exactly the same as the relationship between the forward-backward
  and Virterbi algorithm. Note that in addition to
  <math|<wide|Q|^><around*|(|j,\<tau\>|)>> we also keep track of the choice
  we made to get the maximum value at each node. This allows us to label the
  presence of the mutation in each node, and hence identify the origin and
  loss points.

  <subsubsection|Results and limitations>

  <big-figure|<image|figures/module_4/patient1_snv_loss_cnv.pdf|400pt|||>|<label|fig:loss_cnv>Results
  of the mutation loss model. The panels show copy number profiles from six
  samples from a patient with high grade serous ovarian cancer. Mutations
  predicted to be lost somewhere in the tree are plotted as sticks at the
  bottom of each panel. Highlighted in yellow is an example of a deletion
  event overlapping predicted losses.>

  With the model in hand the question is how can validate the results of the
  model? One ad-hoc approach we can use is to look at the predicted lost
  mutations. If the intuition that copy number changes are deleting mutations
  is valid, we should expect to see evidence of copy number changes in the
  samples where the mutation is lost. Figure <reference|fig:loss_cnv>
  illustrates this point nicely. Here we see two samples, ROv2 and ROv3,
  where a focal deletion of one allele corresponds to the loss of mutations.

  <big-figure|<image|figures/module_4/phylo1.pdf|285pt|234pt||>|<label|fig:phylo1>Simplified
  example of sample tree construction. On the top we show the phylogeny and
  evolutionary history of mutations. The observed presence/absence data is
  shown on the bottom.>

  <big-figure|<image|figures/module_4/phylo2.pdf|285pt|234pt||>|<label|fig:phylo2>Simplified
  example of sample tree construction. On the top we show the phylogeny and
  evolutionary history of mutations. The observed presence/absence data is
  shown on the bottom.>

  <big-figure|<image|figures/module_4/phylo3.pdf|285pt|234pt||>|<label|fig:phylo3>Simplified
  example of sample tree construction. On the top we show the phylogeny and
  evolutionary history of mutations. The observed presence/absence data is
  shown on the bottom.>

  No model is perfect, and the model described so far has one critical
  weakness. If any sample is balanced mixture of two clones then it will fail
  badly. To understand the issue we need to consider a hidden assumption we
  are making. Namely, the probability of clonal presence is accurate. If we
  assign a high probability of presence to sub-clonal mutations from clones
  from different parts of the tree we have a problem.\ 

  To understand this issue first consider Figure <reference|fig:phylo1> which
  illustrates an ideal case of building sample trees. There is no mutation
  loss in this example. In addition the no samples are mixtures of clones
  from different parts of the tree. Thus we can look at the observed
  presence/absence data and easily determine samples S1 and S2 are more
  similar to each other than S3. Now consider Figure <reference|fig:phylo2>
  where sample S2 is a mixture of clones from different parts of the tree.
  The presence/absence data is know ambigous since S2 shares an equal number
  of mutations with S1 and S3. Finally consider Figure <reference|fig:phylo3>
  where we have mutation loss in sample S2. Now sample S1 shares the same
  number of mutations with S2 as S3. The model we constructed accounts for
  the case in Figure <reference|fig:phylo3> but not <reference|fig:phylo2>.
  As result, when we do have mixtures of divergent lineages in a sample the
  model is forced to use mutation loss to explain the mutations which do not
  fit the tree. This can lead to major failures for the method as illustrated
  in Figure <reference|fig:loss_error>. In this example we see a large number
  mutations predicted to be lost. Morevover, they are uniformly spread across
  the genome with no corroborating copy number changes is most cases. The
  problem is that the second sample, LOv2, is actually a mixture of clones
  from LOv1 and the other samples. <big-figure|<image|figures/module_4/patient_9_snv_loss_cnv.pdf|400pt|||>|<label|fig:loss_error>Example
  of erroneous model predictions.>

  The main point of this discussion is to illustrate that no model is
  perfect. It is important to understand the limitations of any model, and
  critically inspect the results. In the study where we first used this
  model, we were well aware of the deficiencies of the model. To address them
  we acknowledged we could not apply this approach to all patients in the
  study. For the cases we could not, we reverted to using single cell
  sequencing which was significantly more expensive and time consuming.

  <subsection|Discussion>

  In this module we discussed proabilistic phylogenetic models. We looked at
  the data types that are available and discussed the challenges associated
  with them. We then reviewed the basics of probabilistic phylogenetic
  models. Finally we developed a novel model for inferring phylogenies from
  bulk sequencing. We saw the success and limitations of this model when
  applied to real data.

  \;
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-type|letter>
    <associate|page-width|auto>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-10|<tuple|1.3.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-11|<tuple|1.3.4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-12|<tuple|1.4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-13|<tuple|1.4.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-14|<tuple|1.4.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-15|<tuple|1.4.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-16|<tuple|1.4.4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-17|<tuple|1.4.5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-18|<tuple|1.4.6|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-19|<tuple|1.4.7|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-2|<tuple|1.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-20|<tuple|1.4.8|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-21|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-22|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-23|<tuple|3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-24|<tuple|4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-25|<tuple|5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-26|<tuple|1.5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-3|<tuple|1.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-4|<tuple|1.2.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-5|<tuple|1.2.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-6|<tuple|1.2.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-7|<tuple|1.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-8|<tuple|1.3.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-9|<tuple|1.3.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:loss_cnv|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:loss_error|<tuple|5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:phylo1|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:phylo2|<tuple|3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:phylo3|<tuple|4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<surround|<hidden|<tuple>>||Results of the mutation loss
      model. The panels show copy number profiles from six samples from a
      patient with high grade serous ovarian cancer. Mutations predicted to
      be lost somewhere in the tree are plotted as sticks at the bottom of
      each panel. Highlighted in yellow is an example of a deletion event
      overlapping predicted losses.>|<pageref|auto-21>>

      <tuple|normal|<surround|<hidden|<tuple>>||Simplified example of sample
      tree construction. On the top we show the phylogeny and evolutionary
      history of mutations. The observed presence/absence data is shown on
      the bottom.>|<pageref|auto-22>>

      <tuple|normal|<surround|<hidden|<tuple>>||Simplified example of sample
      tree construction. On the top we show the phylogeny and evolutionary
      history of mutations. The observed presence/absence data is shown on
      the bottom.>|<pageref|auto-23>>

      <tuple|normal|<surround|<hidden|<tuple>>||Simplified example of sample
      tree construction. On the top we show the phylogeny and evolutionary
      history of mutations. The observed presence/absence data is shown on
      the bottom.>|<pageref|auto-24>>

      <tuple|normal|<surround|<hidden|<tuple>>||Example of erroneous model
      predictions.>|<pageref|auto-25>>
    </associate>
    <\associate|toc>
      1.<space|2spc>Phylogenetic analysis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Overview
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2.<space|2spc>Data types
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|1.2.1.<space|2spc>Bulk sequencing
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.2.2.<space|2spc>Single cell sequencing
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.2.3.<space|2spc>Summary
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|1.3.<space|2spc>Probabilistic phylogenetic
      models <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|2tab>|1.3.1.<space|2spc>Tree topology
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|1.3.2.<space|2spc>Transition probabilities
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|2tab>|1.3.3.<space|2spc>Tree probability
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|2tab>|1.3.4.<space|2spc>Bayesian phylogenetics
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|1tab>|1.4.<space|2spc>Probabilistic model for
      mutation loss <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|1.4.1.<space|2spc>Problem statement
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|2tab>|1.4.2.<space|2spc>Model description
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|2tab>|1.4.3.<space|2spc>Probability of clonal
      presence <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|2tab>|1.4.4.<space|2spc>Tree notation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|2tab>|1.4.5.<space|2spc>Tree probability
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|2tab>|1.4.6.<space|2spc>Inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|2tab>|1.4.7.<space|2spc>Inferring origin,
      presence and loss of mutations <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|2tab>|1.4.8.<space|2spc>Results and limitations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>
    </associate>
  </collection>
</auxiliary>