<TeXmacs|1.99.8>

<style|<tuple|tmbook|british>>

<\body>
  <section|Copy number variation>

  In the last module we looked at probabilistic models for SNV data. In this
  module we will switch gears and consider copy number variants (CNVs). The
  main difference between the two data types from a modelling perspective is
  spatial correlation. If we look at total read counts between adjacent
  positions in the genome they will be strongly correlated. This requires us
  to consider models which can account for the spatial correlation. In
  contrast, it fairly reasonable to model SNVs as independent observations.

  In this module we will explore how to model CNVs. We start by discussing
  some data representations. Next we will construct a simple model for
  inferring total copy number (TCN) from homogeneous cell populations. We
  will then look at some of the practical limitations of this model and how
  they can be addressed.\ 

  <subsection|Data representation>

  Again we will assume we have aligned data from a high throughput sequencing
  experiment. When dealing with CNVs there are a few different ways to
  summarise the aligned read data for modelling purposes. The obvious one is
  to use allele counts in the same way as SNV data. The main downside of this
  approach is that we will typically be interested in every position in the
  genome, whereas we usually consider far fewer positions when dealing with
  SNVs. This can be computationally expensive, and make model fitting
  impractically slow. One solution is to restrict the set of positions we
  consider. A common choice is to look at heterozygous SNPs identified from
  matched normal tissue. This reduces the number of positions from billions
  to millions, which is significantly more tractable. The other benefit of
  heterozygous SNPs is that they can be used infer allele specific copy
  number as we will discuss later.

  An alternative way to represent the data is to divide the genome into bins,
  typically of equal size though this is not essential. We can then count the
  total number of reads within the bin and use this as our data. This
  approach is attractive computationally as we typically use bins with
  lengths on the order of <math|10<rsup|3>-10<rsup|6>> which means we will
  have between <math|10<rsup|3>-10<rsup|6>> data points. Beyond the
  computational savings, binning the data can also be useful when the
  sequencing coverage is low. This has become a more relevant consideration
  in recent years with the development of single cell whole genome sequencing
  platforms. Single cells are often sequenced at low coverage (typically
  0.001x-0.1x) as opposed to the higher coverage of bulk (typically
  30x-100x). The main reason to do this is cost, though technical issues such
  as library diversity also factor in. The most delicate issue in this case
  is choosing an appropriate bin size.

  A variation on the binning strategy is to use bins defined by haplotype
  blocks. A haplotype block is segment of the genome where we are able to
  <with|font-shape|italic|phase> heterozygous SNPs. The point of using these
  blocks is to combine the read counts from a larger set of heterozygous SNPs
  in order to improve allele specific copy number inference. For each block
  we then report three values: the total read count, the number of reads
  supporting the <math|a> allele at heterozygous positions, and the number of
  reads supporting the <math|b> allele at heterozygous positions.

  <subsection|Modelling spatial correlation>

  Regardless of the data representation, we need to consider spatial
  correlation when modelling CNV data. The rational is that copy number
  changes typically affect large regions of the genome that will span many
  data points. Thus it more likely the copy number of adjacent points is the
  same. There are a few commonly used strategies to deal with this problem.

  The simplest strategy is to <with|font-shape|italic|segment>, that is
  identify regions in the same state, before performing any other analysis.
  This can be readily done by standard tools such as circular binary
  segmentation (CBS). Once the data has been segmented we can ignore the
  spatial correlation and treat the data as independent. This is the approach
  used by tools such as ASCAT and the battenberg algorithm. This makes the
  modelling task easier and generally leads to faster methods. The downside
  is that we cannot leverage additional information from the model during the
  segmentation step. This in turn means that we cannot improve the
  segmentation as we learn other features of the data such as the tumour
  content.

  The other major approach is to use probabilistic models with spatial
  correlation. By far the most common approach in this category is to use
  hidden Markov models (HMM). Alternative models such as Kalman filters and
  more general state space models are not widely used in cancer genomics. The
  appeal of HMMs is that they are relatively simple but still able to capture
  some notion of spatial correlation. The key assumption is that a data point
  only depends on the point immediately preceding it. This is the Markov
  assumption of the model and is important for developing efficient inference
  algorithms. The main downside of HMMs is that they have a fairly limited
  ability to model <with|font-shape|italic|state duration>. Informally this
  means that very long segments tend to be discouraged. This will manifest as
  the inference of rapid changes from a state and back into it. Despite this,
  HMMs generally work well and will be the focus of the remainder of this
  module.

  <subsection|Hidden Markov models>

  <subsubsection|Description>

  An HMM is a latent variable model, where the latent (or hidden) variables
  have a dependency structure. The latent variables take values in some
  discrete state space <math|\<cal-X\>>. The observed values depend on the
  latent variables in much the same way as a mixture model. In fact an
  alternative way to define an HMM is as a dynamic mixture model, where the
  probability of belonging to a cluster changes as we move along the
  positions. To define an HMM we need three pieces:

  <\enumerate-numeric>
    <item><math|\<b-pi\>> - The initial state vector. This is a vector of
    length <math|<around*|\||\<cal-X\>|\|>> where <math|\<pi\><rsub|x>> is
    the probability the first hidden variable in the sequence takes value
    <math|x>.

    <item><math|A> - The transition matrix. This is a
    <math|<around*|\||\<cal-X\>|\|>\<times\><around*|\||\<cal-X\>|\|>> matrix
    where <math|A<rsub|i j>> is the probability that the current hidden value
    takes on state <math|j> given the previous one was in state <math|i>.

    <item><math|F> - The emission distribution. This is the distribution for
    the observed data which depends on the associated hidden state. The
    hidden state is thus like the cluster indicator in mixture models,
    selecting which paramter is used for <math|F> to generate the data point.
  </enumerate-numeric>

  Here we assume that <math|K=<around*|\||\<cal-X\>|\|>> is known and fixed.
  This assumption can be relaxed in much the same way as mixture models using
  non-parametric Bayesian priors. However, the resulting models are usually
  computationally demanding to fit.

  <\remark>
    \ Though non-parametric HMMs have not been widely used for studying
    cancer, they could solve some issues. First, like we did in SNVs we could
    use them to infer the number of clones. Second, we could relax the
    maximum copy number constraint that HMM based methods require.
  </remark>

  The basic model for a Bayesian HMM is as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-pi\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-gamma\><rsub|k>|)>>>|<row|<cell|z<rsub|1>\|\<b-pi\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-pi\>|)>>>|<row|<cell|z<rsub|n>\|z<rsub|n-1>,A>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|A<rsub|z<rsub|n-1>\<cdot\>>|)>>>|<row|<cell|\<theta\><rsub|k>>|<cell|\<sim\>>|<cell|G>>|<row|<cell|x<rsub|n>\|z<rsub|n>,\<b-theta\>>|<cell|\<sim\>>|<cell|F<around*|(|\<cdot\>\|\<theta\><rsub|z<rsub|n>>|)>>>>>
  </eqnarray>

  where <math|A<rsub|k\<cdot\>>> is used to denote the <math|k<rsup|th>> row
  of <math|A>, <math|F> and <math|G> are arbitrary distributions with
  densities <math|f> and <math|g>. We can write down the joint distribution
  for this model as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>,\<b-z\>,\<b-pi\>,A,\<b-theta\>|)>>|<cell|=>|<cell|<wide*|p<around*|(|\<b-pi\>|)><big|prod><rsub|k=1><rsup|K>p<around*|(|\<theta\><rsub|k>|)>
    p<around*|(|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>|)>|\<wide-underbrace\>><rsub|prior>\<times\><wide*|p<around*|(|x<rsub|1>\|z<rsub|1>,\<b-theta\>|)>p<around*|(|z<rsub|1>|)><big|prod><rsub|n=2><rsup|N>p<around*|(|x<rsub|n>\|z<rsub|n>,\<b-theta\>|)>p<around*|(|z<rsub|n>\|z<rsub|n-1>|)>|\<wide-underbrace\>>>>|<row|<cell|>|<cell|=>|<cell|p<around*|(|\<b-pi\>|)><big|prod><rsub|k=1><rsup|K>p<around*|(|\<theta\><rsub|k>|)>
    p<around*|(|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>|)>\<times\>>>|<row|<cell|>|<cell|>|<cell|<big|prod><rsub|k=1><rsup|K><around*|[|\<pi\><rsub|k>
    f<around*|(|x<rsub|1>\|\<theta\><rsub|k>|)>|]><rsup|\<bbb-I\><around*|(|z<rsub|1>=k|)>><big|prod><rsub|n=2><rsup|N><big|prod><rsub|\<ell\>=1><rsup|K>f<around*|(|x<rsub|n>\|\<theta\><rsub|\<ell\>>|)><rsup|\<bbb-I\><around*|(|z<rsub|1>=k|)>><big|prod><rsub|k=1><rsup|K>A<rsub|k
    \<ell\>><rsup|\<bbb-I\><around*|(|z<rsub|n-1>=k,z<rsub|n>=\<ell\>|)>
    >>>>>
  </eqnarray>

  <subsubsection|Inference>

  One of the most useful features of HMMs is that we can efficiently compute
  the conditional distribution of the hidden states. This is done using the
  forward-backward (FB) algorithm, which is a dynamic programming algorithm.
  Applying the FB we can easily implement an expectation maximisation
  algorithm to estimate the MAP parameters. Once we have those parameters we
  can use another recursive algorithm, the Viterbi algorithm, to find the
  most probable sequence of hidden states.

  <\remark>
    The forward and backward recursions are also useful when using MCMC
    methods. They can be used to develop a simple Gibbs sampler that
    sequentially updates the hidden states. Another approach is to use the
    <with|font-shape|italic|forward filtering backward sampling> algorithm
    which updates the entire chain.
  </remark>

  In what follows we will suppress the dependencies on the model parameters
  to keep the notation simple. We will just state the required recursion
  formulas here. Interested readers can find the detailed derivation in [cite
  bishop]. We have two recursions the forward,
  <math|\<alpha\><around*|(|z<rsub|n>|)>>, and the reverse,
  (<math|\<beta\><around*|(|z<rsub|n>|)>>).

  <\eqnarray>
    <tformat|<table|<row|<cell|\<alpha\><around*|(|z<rsub|n>=k|)>>|<cell|=>|<cell|p<around*|(|x<rsub|n>\|z<rsub|n>=k|)>
    <big|sum><rsub|\<ell\>=1><rsup|K>\<alpha\><around*|(|z<rsub|n-1>=l|)>A<rsub|\<ell\>
    k> >>|<row|<cell|\<beta\><around*|(|z<rsub|n>=k|)>>|<cell|=>|<cell|<big|sum><rsub|\<ell\>=1><rsup|K>\<beta\><around*|(|z<rsub|n+1>=\<ell\>|)>
    p<around*|(|x<rsub|n+1>\|z<rsub|n+1>=\<ell\>|)> A<rsub|k \<ell\>> >>>>
  </eqnarray>

  The main use of these recursions is to compute the following quantities

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>|)>>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|K>\<alpha\><around*|(|z<rsub|N>=k|)>>>|<row|<cell|\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n>=k|)>|]>>|<cell|=>|<cell|<frac|p<around*|(|z<rsub|n>=k\|\<b-x\>|)>|p<around*|(|\<b-x\>|)>>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<alpha\><around*|(|z<rsub|n>=k|)>\<beta\><around*|(|z<rsub|n>=k|)>|p<around*|(|\<b-x\>|)>>>>|<row|<cell|\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n-1>=k,z<rsub|n>=\<ell\>|)>|]>>|<cell|=>|<cell|p<around*|(|z<rsub|n-1>=k,z<rsub|n>=\<ell\>\|\<b-x\>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<alpha\><around*|(|z<rsub|n-1>-k|)>p<around*|(|x<rsub|n-1>\|z<rsub|n-1>=k|)>p<around*|(|x<rsub|n>\|z<rsub|n>=\<ell\>|)>\<beta\><around*|(|z<rsub|n>=\<ell\>|)>|p<around*|(|\<b-x\>|)>>>>>>
  </eqnarray>

  These become useful in the EM algorithm during the maximisation step when
  we seek to optimise <math|log p<around*|(|\<b-x\>,\<b-pi\>,A,\<b-theta\>|)>>.

  <\eqnarray>
    <tformat|<table|<row|<cell|log p<around*|(|\<b-x\>,\<b-pi\>,A,\<b-theta\>|)>>|<cell|=>|<cell|\<bbb-E\><rsub|\<b-z\>><around*|[|log
    p<around*|(|\<b-x\>,\<b-z\>,\<b-pi\>,A,\<b-theta\>|)>|]>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|K>\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|1>=k|)>|]>log
    \<pi\><rsub|k>+>>|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|n=1><rsup|N><big|sum><rsub|k=1><rsup|K>\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n>=k|)>|]>log
    f<around*|(|x<rsub|n>\|\<theta\><rsub|k>|)>+>>|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|n=2><rsup|N><big|sum><rsub|k=1><rsup|K><big|sum><rsub|\<ell\>=1><rsup|K>\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n-1>=k,z<rsub|n>=\<ell\>|)>|]>log
    A<rsub|k \<ell\>>+>>|<row|<cell|>|<cell|>|<cell|log
    p<around*|(|\<b-pi\>\|\<b-kappa\>|)>+<big|sum><rsub|k=1><rsup|K>log
    p<around*|(|A<rsub|k\<cdot\>>\|\<b-gamma\>|)>+<big|sum><rsub|k=1><rsup|K>log
    p<around*|(|\<theta\><rsub|k>|)>>>>>
  </eqnarray>

  This optimisation takes a closed form for <math|A> and <math|\<b-pi\>>

  <\eqnarray>
    <tformat|<table|<row|<cell|<wide|\<pi\>|^><rsub|k>>|<cell|\<propto\>>|<cell|\<kappa\><rsub|k>+\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|1>=k|)>|]>>>|<row|<cell|<wide|A|^><rsub|k
    \<ell\>>>|<cell|\<propto\>>|<cell|\<gamma\><rsub|k
    \<ell\>>+<big|sum><rsub|k=1><rsup|K><big|sum><rsub|\<ell\>=1><rsup|K>\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n-1>=k,z<rsub|n>=\<ell\>|)>|]>>>>>
  </eqnarray>

  Computing the MAP value of the parameters <math|\<theta\><rsub|k>> will not
  be a closed form in general so often numerically optimise using the
  gradient given by

  <\eqnarray>
    <tformat|<table|<row|<cell|<frac|\<partial\>|\<partial\>
    \<theta\><rsub|k>>log p<around*|(|\<b-x\>,\<b-pi\>,A,\<b-theta\>|)>>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N>\<bbb-E\><rsub|\<b-z\>><around*|[|\<bbb-I\><around*|(|z<rsub|n>=k|)>|]><frac|\<partial\>|\<partial\>\<theta\><rsub|k>>log
    f<around*|(|x<rsub|n>\|\<theta\><rsub|k>|)>+<frac|\<partial\>|\<partial\>\<theta\><rsub|k>>log
    p<around*|(|\<theta\><rsub|k>|)>>>>>
  </eqnarray>

  <\remark>
    Thus far we consider only a single sequence. We will see in the next
    section that we will treat each chromosome as a separate sequence. The
    above equations are simple to modify to accommodate this variation.
  </remark>

  <subsection|Homogeneous sample model>

  We now turn our attention to model building. Here we will consider the
  problem of inferring copy number profiles from a sample that has no clonal
  population structure and no normal contamination. This model would also
  apply equally well to single cell data. We assume that the data is binned
  total read counts. Given this we will not try to infer allele specific copy
  number, but instead focus on total copy number. Let's introduce some
  notation. Let\ 

  <\itemize-dot>
    <item><math|m> index the chromosomes

    <item><math|M> be the number of chromosomes

    <item><math|n> index the bins in a chromosome

    <item><math|N<rsub|m>> be the number of bins for chromosome <math|m>

    <item><math|x<rsup|m><rsub|n>> denote the number of reads in the
    <math|n<rsup|th>> bin on chromosome <math|m>

    <item><math|z<rsub|n><rsup|m>> denote the hidden state for the
    <math|n<rsup|th>> bin on chromosome <math|m>
  </itemize-dot>

  <subsubsection|Model description>

  The first step is think about an appropriate distribution to model the read
  counts. Given the observed data takes non-negative integer values, a simple
  choice is the Poisson distribution. Next we need to think about our hidden
  states. We would like them to encode the copy number of a bin, so we will
  assume they take values in the set <math|<around*|{|1,\<ldots\>,C|}>> where
  <math|C> is some maximum possible copy number.

  <\remark>
    We will ignore homozygous deletions for the simplicity.
  </remark>

  Now we need to think about the prior distribution for the parameters of the
  Poisson data distribution. The parameter of a Poisson needs to be a
  positive real value. A very common choice for this type of parameter is the
  Gamma distribution. If we sample the parameter for each state from an
  arbitrary Gamma there is no relationship between the copy number and data.
  We will take a different approach here. We will sample a single value
  <math|r> from a Gamma and let the Poisson parameter for state <math|c> be
  <math|\<theta\><rsub|c>=c r>. The parameter <math|r> can be interpreted as
  the haploid bin coverage i.e. the coverage of a bin with copy number 1.
  Thus we are assuming the total number of reads in a bin is a linear
  function of copy number. We can set the parameters <math|a> and <math|b> in
  the Gamma prior so the distribution has a mean equal to the haploid depth
  we expect. Alternatively they can be set to values which lead to a
  <with|font-shape|italic|vague> distribution i.e. one with high variance.

  We finish off by setting the prior parameters for <math|\<b-pi\>> and
  <math|A>. There is no obvious reason to prefer any particular initial copy
  number on a chromosome, so we let <math|\<b-kappa\>=<around*|(|1,\<ldots\>,1|)>>
  which leads to the uniform prior on the simplex (vectors that sum to one).
  We would like to encode a preference for adjacent bins to share the same
  copy number. To achieve this we let <math|\<b-gamma\><rsub|k
  \<ell\>>=\<zeta\>> for <math|\<ell\>\<neq\>k> and <math|\<gamma\><rsub|k
  k>=\<zeta\>+\<xi\>> for some values <math|\<zeta\>> and <math|\<xi\>>. To
  be concrete we will let <math|\<zeta\>=1> and <math|\<xi\>=10>. This prior
  will favour the hidden variables to stay in the same state as their
  neighbour. The full model is then

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-pi\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-gamma\><rsub|k>|)>>>|<row|<cell|z<rsub|1><rsup|m>\|\<b-pi\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-pi\>|)>>>|<row|<cell|z<rsub|n><rsup|m>\|z<rsub|n-1><rsup|m>,A>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|A<rsub|z<rsub|n-1><rsup|m>\<cdot\>>|)>>>|<row|<cell|r>|<cell|\<sim\>>|<cell|<text|Gamma><around*|(|\<cdot\>\|a,b|)>>>|<row|<cell|\<theta\><rsub|c>\|r>|<cell|=>|<cell|r
    c>>|<row|<cell|x<rsub|n><rsup|m>\|z<rsub|n><rsup|m>,\<b-theta\>>|<cell|\<sim\>>|<cell|Poisson<around*|(|\<cdot\>\|\<theta\><rsub|z<rsub|n><rsup|m>>|)>>>>>
  </eqnarray>

  <subsubsection|Inference>

  We can use the EM algorithm as discussed above. Minor modification need to
  be made to account for the fact we have multiple sequences (chromosomes).

  <subsubsection|Limitations>

  The main practical issue this model will run into is the assumption that
  the read counts only depend on the total copy number number. In reality
  issues such as mappability and GC content will also impact read depth. One
  approach is to perform some form of normalisation to the read counts. This
  typically changes the values to be non-integer, so we would need to abandon
  the Poisson assumption. A standard way forward is then to work with log
  transformed counts. We can then normalise those and model the data as
  following a Normal distribution instead of a Poisson.

  If we do not want to abandon the current model, then we could try
  incorporating the additional data about covariates such as GC content into
  the model. Let <math|g<rsub|n><rsup|m>> denote some measure of the GC
  content for bin <math|n> on chromosome <math|n>. A very simple way to
  incorporate this information is to modify the data distribution as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|x<rsub|n><rsup|m>\|g<rsub|n><rsup|m>,z<rsub|n><rsup|m>,\<b-theta\>>|<cell|\<sim\>>|<cell|Poisson<around*|(|\<cdot\>\|g<rsub|n><rsup|m>\<theta\><rsub|z<rsub|n><rsup|m>>|)>>>>>
  </eqnarray>

  We have simply scaled the Poisson parameter by the GC content information.
  Often incorporating observed data directly into the likelihood causes
  problems. One reason is the observed data maybe on the wrong scale or it
  may be error prone itself. A simple solution is to introduce an additional
  layer to the model. For example we could pick some distribution and set its
  mean to <math|g<rsub|n><rsup|m>> denoted by
  <math|H<around*|(|\<cdot\>\|g<rsub|n><rsup|m>|)>>. Then we augment the
  model as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|h<rsub|n><rsup|m>>|<cell|\<sim\>>|<cell|H<around*|(|\<cdot\>\|g<rsub|n><rsup|m>|)>>>|<row|<cell|x<rsub|n><rsup|m>\|h<rsub|n><rsup|m>,z<rsub|n><rsup|m>,\<b-theta\>>|<cell|\<sim\>>|<cell|Poisson<around*|(|\<cdot\>\|h<rsub|n><rsup|m>
    \<theta\><rsub|z<rsub|n><rsup|m>>|)>>>>>
  </eqnarray>

  Now the observed value <math|g<rsub|n><rsup|m>> only impacts the likelihood
  indirectly, and we could tune additional parameters for <math|H> to fit the
  data better. A reasonable choice for <math|H> in this case would again be a
  Gamma distribution with appropriate parameters.

  The other major deficiency of the model is the assumption of a Poisson
  likelihood. The Poisson only has a single parameter which controls its mean
  and variance. As a result it can often be a poor fit to real data that has
  more variability than it can model. This is called overdispersion and we
  touched on in the first module. The standard solution is to switch to an
  overdispersed distribution such as the Negative-Binomial. We keep the mean
  of this distribution the same as for the Poisson, but introduce an
  additional variance parameter. It may be useful to introduce separate
  variance parameters for each copy number state, or parameterise the
  variance to depend on the state in way analogous to the mean.

  <subsection|Non-homogeneous sample model>

  The previous model was quite simplistic in assuming a homogeneous
  population of cells. In patient data the most major violation of this
  assumption would be normal contamination. Clonal population structure will
  also be an issue. In this section we will look at how we can improve the
  basic model to address these issues.

  <subsubsection|Normal contamination>

  Normal contamination is actually relatively simple to address. Let <math|t>
  be the tumour content of the sample. Then for any bin we have proportion
  <math|t> of cells with copy number <math|c> and proportion
  <math|<around*|(|1-t|)>> of cells with copy number 2 (assuming autosomes).
  We do not know the tumour content so we will need to specify a prior
  distribution. A simple choice is of course the continuous Uniform
  distribution. However, an informative Beta distribution maybe more useful
  if there is information from pathology estimates or other sources. Here we
  will use the simple Uniform.

  The updated model is then as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-pi\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-gamma\><rsub|k>|)>>>|<row|<cell|z<rsub|1><rsup|m>\|\<b-pi\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-pi\>|)>>>|<row|<cell|z<rsub|n><rsup|m>\|z<rsub|n-1><rsup|m>,A>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|A<rsub|z<rsub|n-1><rsup|m>\<cdot\>>|)>>>|<row|<cell|r>|<cell|\<sim\>>|<cell|<text|Gamma><around*|(|\<cdot\>\|a,b|)>>>|<row|<cell|t>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>>|<row|<cell|\<theta\><rsub|c>\|r,t>|<cell|=>|<cell|r<around*|(|2
    <around*|(|1-t|)>+t \ c|)>>>|<row|<cell|x<rsub|n><rsup|m>\|z<rsub|n><rsup|m>,\<b-theta\>>|<cell|\<sim\>>|<cell|Poisson<around*|(|\<cdot\>\|\<theta\><rsub|z<rsub|n><rsup|m>>|)>>>>>
  </eqnarray>

  Inference will get slightly more complicated as we need to infer <math|t>.
  But the same basic steps will work.

  <subsubsection|Ploidy and identifiability>

  One issue we have ignored thus far is whether our model is
  <with|font-shape|italic|identifiable>. Roughly speaking a model is
  identifiable if every combination of parameters maps to a unique likelihood
  value in the frequentist setting, or joint probability value in the
  Bayesian setting. If a model is unidentifiable then we can have multiple
  MLE or MAP solutions. This raises the issue of how to select the best
  solution. In a fully Bayesian analysis, unidenitifiability is not
  theoretically an issue. We will simply observe a multi-modal posterior
  distribution. In practice unidenitifiability can make it hard to implement
  efficient MCMC algorithms as we need to be able
  <with|font-shape|italic|hop> between modes to explore the posterior. As a
  general rule of thumb it is best to avoid constructing models which are
  unidentifiable.

  In the case of the simple model we defined for CNV analysis this
  unindentifiability crops up in the term <math|\<theta\><rsub|c>=r c> which
  paramaterises the Poisson observation distribution. The issue we have is
  that we can double the copy number of all segments and half the value of
  <math|r> to obtain the same value of <math|\<theta\><rsub|c>>. The
  likelihood of these solutions is thus identical. Because we have a prior on
  <math|r> the joint probability of these two different solution may be
  slightly different. However, the data is not informing the solution only
  the prior. In this case we say the model is <with|font-shape|italic|weakly
  identifiable>. This situation is worrying because we need to trust the
  prior for <math|r> if we are to trust the MAP solution. If we adopted the
  frequentist viewpoint, things would worse still as we would have no way to
  choose a solution. It may seem that adding tumour content into the model
  can help with this issue as our value of <math|\<theta\>> now depends on
  <math|t> as well. Unfortunately this is not the case, as we can also adjust
  the value of <math|t> to derive equally likely solutions.\ 

  In the field of CNV analysis the term <with|font-shape|italic|ploidy> is
  loosely used to refer to the average copy number of a sample. The exact
  definition varies slightly depending on the model, but the basic idea is
  the same. To the best of the author's knowledge all models for CNV analysis
  have the problem of unidentifiability caused by ploidy i.e. scaling the
  total copy number of all segments. The common approach employed is to
  manually select a solution from the of MLE or MAP estimates that have been
  found. This issue afflicts both bulk and single CNV analysis. Automating
  the selection of ploidy solution remains an open problem, and in the
  absence of any additional information or assumptions on the model appears
  very difficult. This is one reason that many tools for CNV analysis will
  give very different results when analysing the same sample.

  There are a few avenues of research that could be pursued to address this
  issue. If we have additional information about the ploidy of the samples,
  for example from flow cytometry, we could leverage this information. This
  is only applicable to bulk analysis however, and is not full proof. We may
  also have domain specific knowledge about the cancer type we are analysing.
  For example we may know a priori that the genomes are relatively stable, so
  a near diploid solution is preferred. Alternatively we may have cancers
  which are known to undergo genome doubling early, in which case a near
  tetraploid solution would be preferred. Again this is not a general
  solution, as most cancers fall on a spectrum of genomic instability. Some
  single cell sequencing platforms image the cells before sequencing. It is
  possible this could be used to identify cells which are larger and more
  likely to have higher ploidy, though it is unclear whether this will work.
  The final idea would be to model the variance of each copy number state and
  try to extract information from this. Under certain modelling assumptions
  cells with higher overall copy number will have more variability in the
  observed reads counts.

  In summary CNV analysis is challenging because of the ploidy problem.
  Currently there is no full proof automated way to address this issue. In
  practice this means that either manual curation or ad-hoc post-processing
  needs to be performed to identify a solution.

  <subsection|Allele specific copy number model>

  We will briefly outline the steps needed to produce an allele specific copy
  number model. We will assume the data comes in the form of haploid blocks.
  Let <math|x<rsup|a><rsub|n>> be the number of reads supporting the <math|a>
  allele in the block <math|n>, <math|x<rsup|b><rsub|n>> be the number of
  reads supporting the <math|b> allele, <math|x<rsub|n><rsup|t>> be the total
  read depth and <math|x<rsub|n><rsup|l>> the length of the block. Note that
  we no longer assume the blocks are of equal length. This will impact the
  expected number of reads in a block.

  We know let the state space of the HMM be
  <math|<around*|{|<around*|(|c<rsub|a>,c<rsub|b>|)>:c<rsub|a>,c<rsub|b>\<leq\>C|}>>
  i.e. all pairs of allele specific copy numbers. Again we set a maximum
  number of copies <math|C>. Our state space has now grown from size <math|C>
  to <math|C<rsup|2>>, so inference will be much slower.

  We will stick with the Poisson distribution to model total read depth,
  <math|x<rsup|t><rsub|n>>, where we scale the haploid depth by total copy.
  But we will now include the bin length. This alters the interpretation of
  <math|r>, which now is the haploid depth of a position not a bin.

  We also need to model the allele specific counts. A simple choice is to use
  a Binomial for the number reads supporting the <math|b> allele. We will let
  the probability of success be <math|<frac|c<rsub|b>|c<rsub|a>+c<rsub|b>>>
  so the probability depends on the relative copy number of the <math|b>
  allele.

  We suppress the index for chromosomes for clarity. The full model is then
  as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-pi\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|A<rsub|k\<cdot\>>\|\<b-gamma\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-gamma\><rsub|k>|)>>>|<row|<cell|z<rsub|1>\|\<b-pi\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-pi\>|)>>>|<row|<cell|z<rsub|n>\|z<rsub|n-1>,A>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|A<rsub|z<rsub|n-1>\<cdot\>>|)>>>|<row|<cell|r>|<cell|\<sim\>>|<cell|<text|Gamma><around*|(|\<cdot\>\|a,b|)>>>|<row|<cell|\<theta\><rsub|n>\|z<rsub|n>=<around*|(|c<rsub|a>,c<rsub|b>|)>,x<rsub|n><rsup|l>>|<cell|=>|<cell|r
    <around*|(|c<rsub|a>+c<rsub|b>|)> x<rsub|n><rsup|l>>>|<row|<cell|x<rsub|n><rsup|t>\|z<rsub|n>,\<b-theta\>>|<cell|\<sim\>>|<cell|Poisson<around*|(|\<cdot\>\|\<theta\><rsub|n>|)>>>|<row|<cell|x<rsub|n><rsup|b>\|x<rsub|n><rsup|a>,z<rsub|n>=<around*|(|c<rsub|a>,c<rsub|b>|)>>|<cell|\<sim\>>|<cell|<text|Binomial><around*|(|\<cdot\>\|x<rsub|n><rsup|a>+x<rsub|n><rsup|b>,<frac|c<rsub|b>|c<rsub|a>+c<rsub|b>>|)>>>>>
  </eqnarray>

  Fitting this model is no harder than the original model, as we have not
  gained any additional parameters. Evaluating the likelihood will become
  more expensive due to the added Binomial term.

  <subsection|Discussion>

  <subsubsection|Summary>

  In this module we saw how to construct probabilistic models for CNV data.
  We discussed different data representations that are commonly used. We then
  explored using HMM models to capture the spatial correlation of the data.
  We constructed a simple model for inferring total copy number profiles from
  homogeneous populations. We then explored how this model could be altered
  to address more complex issues. We note that all of these approaches could
  be combined to produce a more general model.

  <subsubsection|Other approaches>

  We focused on using HMMs in this module, primarily to illustrate how to
  handle spatial correlation in data. We could have avoided the use of HMMs,
  and the computational complexity this entails, by segmenting the data
  first. This is a widely used approach the effectively makes the data
  independent. This data can then easily be modelled using ideas similar to
  those in the previous module.

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
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|1.4.3|5>>
    <associate|auto-11|<tuple|1.5|5>>
    <associate|auto-12|<tuple|1.5.1|5>>
    <associate|auto-13|<tuple|1.5.2|6>>
    <associate|auto-14|<tuple|1.6|6>>
    <associate|auto-15|<tuple|1.7|6>>
    <associate|auto-16|<tuple|1.7.1|6>>
    <associate|auto-17|<tuple|1.7.2|?>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|1>>
    <associate|auto-4|<tuple|1.3|2>>
    <associate|auto-5|<tuple|1.3.1|2>>
    <associate|auto-6|<tuple|1.3.2|2>>
    <associate|auto-7|<tuple|1.4|4>>
    <associate|auto-8|<tuple|1.4.1|4>>
    <associate|auto-9|<tuple|1.4.2|4>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      1.<space|2spc>Copy number variation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Data representation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2.<space|2spc>Modelling spatial
      correlation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|1.3.<space|2spc>Hidden Markov models
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.3.1.<space|2spc>Description
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.3.2.<space|2spc>Inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|1.4.<space|2spc>Homogeneous sample model
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|2tab>|1.4.1.<space|2spc>Model description
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|1.4.2.<space|2spc>Inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|2tab>|1.4.3.<space|2spc>Limitations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|1.5.<space|2spc>Non-homogeneous sample
      model <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|1.5.1.<space|2spc>Normal contamination
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|1.5.2.<space|2spc>Ploidy and
      identifiability <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|1tab>|1.6.<space|2spc>Allele specific copy number
      model <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|1.7.<space|2spc>Discussion
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|2tab>|1.7.1.<space|2spc>Summary
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|2tab>|1.7.2.<space|2spc>Other approaches
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>
    </associate>
  </collection>
</auxiliary>