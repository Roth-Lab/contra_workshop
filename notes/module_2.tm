<TeXmacs|1.99.13>

<style|<tuple|tmbook|british|old-dots>>

<\body>
  <section|Inferring clonal population structure from SNV data>

  In this module we consider the problem of inferring clonal population
  structure using bulk sequencing data. We specifically consider the case of
  using SNVs as our clonal markers.

  <subsection|Problem>

  We consider the problems of how to use high throughput sequencing to infer
  the clonal population structure of a tumour. This problem is somewhat old
  now, but still remains relevant when looking at large cancer cohort
  datasets that are continuing to be generated. The basic ideas of
  deconvolution are also appearing in other areas so this module should
  provide some useful insight.

  Precisely we will use read count data from SNVs to

  <\enumerate-numeric>
    <item>Infer what proportion of cancer cells harbouring a mutation

    <item>Infer what mutations share the same evolutionary history i.e.
    originate at the same time and are lost in the same subsets of clones
  </enumerate-numeric>

  <big-figure|<image|../figures/module_2/bulk_sequencing_v2.pdf||400pt||>|<label|fig:bulk_sequencing>Schematic
  of bulk sequencing for a tumour experiment. At the top we have the input
  cell population, where stars indicate mutations. In the middle we show the
  aligned reads obtained from performing bulk sequencing. Positions in reads
  are colour coded to match the mutations at the top. Note the proportion of
  reads with a variant is roughly similar to the proportion of input cells
  with the mutation. At the bottom we show the summarised read counts which
  we will use for modelling.>

  We will assume that we have collected one or more tumour samples from a
  patient and performed high throughput sequencing. We assume that SNVs have
  already been identified using a standard variant caller tool. We summarise
  our aligned read count format in terms of three quantities for each SNV
  (see Figure <reference|fig:bulk_sequencing>).

  <\itemize-dot>
    <item><math|a> - The number of reads which match the reference allele.

    <item><math|b> - The number of reads which match the variant allele.

    <item><math|d=a+b> - The total depth of coverage at the locus.
  </itemize-dot>

  To achieve our goals will we need to develop two parts to the model. The
  first part will be a way to account for the effect of <em|mutational
  genotype> and <em|normal contamination>. Mutational genotype refers to the
  fact that not all mutation will be heterozygous diploid events due to
  coincident copy number variation. Normal contamination refers to the fact
  that we also sequence non-malignant cells in real tumours. The second part
  will be a mechanism to cluster the mutations. Here we will use the
  formalism of mixture models. There is one challenging issue, which is that
  we do not know how many clones there are in the sample(s). We will address
  this problem using the Dirichlet process.

  <subsection|Modelling mutation genotype>

  The first issue we tackle is the problem of mutational genotype. Our goal
  in this section is to define a probabilistic model that links the observed
  read count data with cellular prevalence. We can then apply the standard
  Bayesian machinery to compute estimates of cellular prevalence.

  <subsubsection|Background>

  We will make the assumption that a point mutation only occurs once at a
  locus during the evolutionary history of tumour. This is often referred to
  as the <with|font-shape|italic|infinite sites assumption>. This assumption
  is motivated by the fact that the mutation rate is usually relatively low
  compared to the total size of the genome. It can break down at mutational
  hot spots, and this has been observed, for example when multiple
  substitutions are observed at a single locus. So this assumption represents
  our first approximation to truth.

  <big-figure|<image|../figures/module_2/allelic_vs_cellular_frequency.v2.pdf|400pt|||>|<label|fig:ccf_vs_vaf>Example
  of a heterozygous diploid mutation showing why the variant allele frequency
  (VAF) is not the same proportion of cells harbouring the mutation (cellular
  prevalence). On the left we have the input population of cells which all
  have the mutation. On the right we have the observed sequence data where
  only half the reads (on average) have the mutation.>

  The issue we face is that tumour cells are not haploid. In the absence of
  copy number events, they are at least diploid at any autosomal locus.
  Figure <reference|fig:ccf_vs_vaf> illustrates the basic problem. In this
  example we have a set of cells which are all diploid and heterozygous for a
  mutation. All the cells have the mutation, but only half the reads are
  expected to show the mutation. This suggests that if we simply double the
  observed variant allele frequency (VAF) we can get a reasonable estimate of
  the cellular prevalence (proportion of cells with the mutation). The
  problem is of course that cancers are aneuploid. Figure
  <reference|fig:mut_genotype> illustrates the challenges that copy number
  variation poses.

  <big-figure|<image|../figures/module_2/genotype.pdf|400pt|||>|<label|fig:mut_genotype>Effect
  of mutational genotype on observed VAF.>

  <subsubsection|Population structure>

  To overcome the issue of mutational genotype we need to start making a
  model of how the observed data is generated. We break the process down into
  the following steps

  <\enumerate-numeric>
    <item>We select a cell proportional to how prevalent it is in the input
    sample.

    <item>Given the cells genotype we then select a chromosome at random.
    Thus the probability of selecting a mutation is proportional to how many
    copies of the chromosome have the mutation.

    <item>Finally we introduce a small error probability
    <math|\<varepsilon\>>. We assume the probability is the same whether we
    truly sampled the reference and observed an erroneous variant, or vice
    versa. Because of this symmetry the error term only appears in cases when
    all of the chromosomes have the reference or the variant allele.
  </enumerate-numeric>

  <big-figure|<image|../figures/module_2/population_structure.svg|400pt|||>|<label|fig:pop_structure>Illustration
  of the assumed population structure. Here all populations are defined with
  respect to a single mutation. The circular cells are non-malignant and the
  irregularly shaped ones are malignant.>

  Figure <reference|fig:pop_structure> illustrates the basic model. We
  decompose the population into three parts

  <\itemize-dot>
    <item>Normal population - The non-malignant cells

    <item>Reference population - The malignant cells without the mutation

    <item>Variant population - The malignant cells with the mutation
  </itemize-dot>

  These populations are all defined with respect to a single mutation. If we
  were to pick a different mutation the cells which belong to the reference
  and variant population would change. This is a common point of confusion as
  people do not understand how a malignant cell does not have mutations. They
  do have mutations, just not the particular mutation we are thinking about
  at the moment.

  <big-figure|<image|../figures/module_2/population_structure_simple.svg|400pt|||>|<label|fig:pop_structure_simple>Illustration
  of the simplified population structure. In contrast to Figure
  <reference|fig:pop_structure_simple> the mutational genotypes of all cells
  within the reference and variant populations are the same. Note the
  genotypes are different between the populations, as they must be since one
  population has the mutation and the other does not.>

  The point of the population decomposition is that the quantity of interest,
  the cellular prevalence of a mutation, is the proportion of cancer cells in
  the variant population. We are nearly ready to write down the probabilistic
  model but there is one issue. In Figure <reference|fig:pop_structure> we
  assume the genotypes of the malignant cells can be highly variable. While
  realistic, this makes it very hard to write down the probability of
  sampling a variant allele.\ 

  To avoid this issue we make a simplifying assumption which is that the
  mutational genotype of all cells within the variant population is the same.
  We make the same assumption for the reference population as well. Figure
  <reference|fig:pop_structure_simple> illustrates this assumption. We
  emphasise this is an assumption that is likely violated in many cases.
  Ideally we would not need it, but the resulting model would be very complex
  and statistically <with|font-shape|italic|unidentifiable>. One way we will
  try to mitigate the impact of this assumption later when introduce multiple
  samples, is to allow the mutational genotypes to differ between samples.
  Roughly speaking this is then saying we only care about the mutational
  genotype of the dominant clone in the sample.

  <subsubsection|Modelling counts>

  With our simplified decomposition of populations and assumption we can now
  write down the probability of observing a variant allele. To do this we
  need some notation. Let

  <\itemize-dot>
    <item><math|\<phi\>> be the cellular prevalence of the mutation

    <item><math|t> be the proportion of cancer cells in the sample (tumour
    content)

    <item><math|g<rsub|N>> be the genotype of the normal population

    <item><math|g<rsub|R>> be the genotype of the reference population

    <item><math|g<rsub|V>> be the population of the variant population

    <item><math|\<b-psi\>=<around*|(|g<rsub|N>,g<rsub|R>,g<rsub|V>|)>> be the
    vector of genotypes for notational convenience

    <item><math|c<around*|(|g|)>> be the copy number of genotype <math|g>

    <item><math|\<mu\><around*|(|g|)>> be the proportion of chromosomes with
    the mutation for genotype <math|g>. We actually make a small modification
    to accommodate sequence error and let
    <math|\<mu\><around*|(|g|)>=\<varepsilon\>> when there are no mutations
    and <math|\<mu\><around*|(|g|)>=1-\<varepsilon\>> when all the
    chromosomes are mutated.\ 
  </itemize-dot>

  Now we imagine we know the value of all these quantities. We then ask what
  is the probability of sampling a read with the mutation? We imagine that we
  have an infinite pool of cells which have all been lysed (cells broken
  apart) so that we have an infinite mixture DNA fragments. Now the
  probability of sampling a DNA fragment of cell depends of two things: first
  how common that type of cell was, second how many copies of the locus the
  cell has. We can easily derive these from the quantities we have defined.
  The probability of sampling a read from the\ 

  <\itemize-dot>
    <item>normal population is proportional to
    <math|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>>

    <item>reference population is proportional to <math|t
    <around*|(|1-\<phi\>|)> c<around*|(|g<rsub|R>|)>>

    <item>variant population is proportional to <math|t \<phi\>
    c<around*|(|g<rsub|V>|)>>
  </itemize-dot>

  <\remark>
    It maybe surprising that we consider how many copies of the locus a cell
    has when determining the probability of sampling a read from that cell.
    To see why this necessary consider the following thought experiment. What
    would happen if one population had an infinite number of copies of the
    locus? Then it would contribute an infinite number of DNA fragments, and
    we would only sample these.
  </remark>

  <\remark>
    One can derive a slightly different model if assume cells are first
    selected before the lysis step. In this case the number of copies of the
    locus a cell has does not impact the probability of selecting a read from
    the cell. The complication is that we then need to model how we select
    which chromosomes from the cell get sequenced. It is simple if we pick
    one at random, but harder if multiple chromosomes can be sequenced.
  </remark>

  Note we need to normalise the probabilities to sum to one. Formally, let
  <math|E<rsub|i>\<in\><around*|{|N,R,V|}>> be a random variable indicating
  that read <math|i> was sample from population <math|e>. Then

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|E<rsub|i>=e|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|<frac|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>|Z>>|<cell|<text|if>>|<cell|e=N>>|<row|<cell|<frac|t
    <around*|(|1-\<phi\>|)> c<around*|(|g<rsub|R>|)>|Z>>|<cell|<text|if>>|<cell|e=R>>|<row|<cell|<frac|t
    \<phi\> c<around*|(|g<rsub|V>|)>|Z>>|<cell|<text|if>>|<cell|e=V>>>>>>>|<row|<cell|Z>|<cell|=>|<cell|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>+t
    <around*|(|1-\<phi\>|)> c<around*|(|g<rsub|R>|)>+t \<phi\>
    c<around*|(|g<rsub|V>|)>>>>>
  </eqnarray>

  Once we know which population a read comes from, then it is straightforward
  to compute the probability that the read has the mutation. This is simply
  given by <math|\<mu\><around*|(|g<rsub|e>|)>> if <math|E=e>. To formalise
  this let <math|F<rsub|i>> be a random variable indicating if read <math|i>
  is has the mutation. With this notation we have that

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|F<rsub|i>=1\|E<rsub|i>=E|)>>|<cell|=>|<cell|\<mu\><around*|(|g<rsub|e>|)>>>>>
  </eqnarray>

  hence

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|F<rsub|i>=1|)>>|<cell|=>|<cell|<big|sum><rsub|a\<in\><around*|{|N,R,V|}>>p<around*|(|F<rsub|i>=1\|A<rsub|i>=a|)>p<around*|(|A<rsub|i>=a|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>|Z>
    \<mu\><around*|(|g<rsub|N>|)>+<frac|t <around*|(|1-\<phi\>|)>
    c<around*|(|g<rsub|R>|)>|Z> \<mu\><around*|(|g<rsub|R>|)>+<frac|t \<phi\>
    c<around*|(|g<rsub|V>|)>|Z> \<mu\><around*|(|g<rsub|V>|)>>>|<row|<cell|>|<cell|\<assign\>>|<cell|\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>>>>>
  </eqnarray>

  This implies that\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|F<rsub|i>\|\<b-psi\>,\<phi\>,t>|<cell|\<sim\>>|<cell|<text|Bernoulli><around*|(|\<cdot\>\|\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)>>>>>
  </eqnarray>

  Of course we observe more than one read, and the total number of reads with
  a variant is <math|B=<big|sum><rsub|i=1><rsup|d>F<rsub|i>>. There is a
  basic result in probability that says the sum of Bernoulli variables
  follows a Binomial distribution, so we have

  <\eqnarray>
    <tformat|<table|<row|<cell|B\|\<b-psi\>,\<phi\>,t,d>|<cell|\<sim\>>|<cell|<text|Binomial><around*|(|\<cdot\>\|d,\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)>>>>>
  </eqnarray>

  This is sufficient to define the likelihood for the basic model. The only
  thing left is to specify a prior for <math|\<phi\>>. We know that
  <math|\<phi\>\<in\><around*|[|0,1|]>> but very little else for an arbitrary
  mutation. Thus we will use a Uniform continuous prior. So our simple
  Bayesian model is then

  <\eqnarray>
    <tformat|<table|<row|<cell|\<phi\>>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>>|<row|<cell|B\|\<b-psi\>,\<phi\>,t,d>|<cell|\<sim\>>|<cell|<text|Binomial><around*|(|\<cdot\>\|d,\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)>>>>>
  </eqnarray>

  We can then write down the joint distribution

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|B=b,\<b-psi\>,\<phi\>,t,d|)>>|<cell|=>|<cell|p<around*|(|B=b\|\<b-psi\>,\<phi\>,t,d|)>
    p<around*|(|\<phi\>|)>>>|<row|<cell|>|<cell|=>|<cell|<binom|d|b>\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)><rsup|b><around*|(|1-\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)><rsup|d-b>
    \<bbb-I\><around*|(|\<phi\>\<in\><around*|[|0,1|]>|)>>>>>
  </eqnarray>

  and applying Bayes' rule we can get the posterior

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<phi\>\|=b,\<b-psi\>,t,d|)>>|<cell|=>|<cell|<frac|p<around*|(|B=b,\<b-psi\>,\<phi\>,t,d|)>|<big|int>p<around*|(|B=b,\<b-psi\>,\<phi\>,t,d|)>
    \<mathd\>\<phi\>>>>>>
  </eqnarray>

  Now the integral in the bottom does not have a closed form, but it is a one
  dimensional integral and so is trivial to compute numerically. Thus we can
  now compute the posterior distribution for the cellular prevalence of a
  mutation. There is still some work to do, in particular we have assumed we
  know the values of the population genotypes <math|\<b-psi\>> and the tumour
  content <math|t>. We can often get a reasonable estimate of <math|t> from
  other sources such as copy number analysis. However, we are unlikely to
  have any idea about <math|\<b-psi\>>. The genotype of the normal population
  is not an issue, but the other two populations are a mystery. In the next
  section we will discuss how to address this issue.

  <subsubsection|Genotype priors>

  In the previous section we saw that given the mutational genotypes of the
  populations we could infer the cellular prevalence of the mutation. In
  practice we do not have this information. In the Bayesian framework this
  issue appears frequently. We often imagine a model where we have more
  information then we can observe in the form of unobserved or
  <with|font-shape|italic|latent> variables. The solution to this problem is
  quite simple, we place a prior distribution over the variables. Once we
  have done this we can either infer them or
  <with|font-shape|italic|marginalise> them. It is generally preferable to
  marginalise the variables if possible as it reduces the dimensionality of
  the problem. This typically leads to faster inference algorithms as the
  space to be explored is smaller.

  The question we face is how to specify a prior for the genotypes
  <math|\<b-psi\>>? One solution is to place a uniform prior over all
  possible genotypes. However, we would face two difficulties. First, our
  posteriors would be highly uninformative. Since we assume nothing about the
  genotypes, essentially any cellular prevalence would be equally likely
  regardless of the data. The second issue is that this is computationally
  infeasible as we would need to sum over an infinite set.

  <big-figure|<image|../figures/module_2/copy_number_with_mutation.png|400pt|||>|<label|fig:cn_snv>Illustration
  of allele specific copy number profile. The red line is the major copy
  number and the blue line is the minor copy number. The star indicates the
  location of an SNV. The major copy number of this SNV is 2 and the minor is
  1.>

  To address this issue we will make use of some auxiliary information that
  we typically have, namely the copy number of the locus. We will assume that
  <with|font-shape|italic|allele specific> copy number profiles are available
  for the samples. These can be generated from micro-array, whole genome
  sequencing or exome sequencing. In the next module we will discuss the
  details of how these are actually generated. Figure <reference|fig:cn_snv>
  provides a schematic example of the information. Here we have an SNV and
  the copy number profile. We can see from this figure that the major copy
  number of the SNV is 2 and the minor copy is 1.

  <\remark>
    We use the terminology major/minor copy number to refer to the parental
    chromosome with more/less copies.
  </remark>

  <\remark>
    Many tools for copy number profiling of cancers provide sub-clonal copy
    number predictions. We ignore this complication, and simply take the
    profile of the most prevalent copy number clone.
  </remark>

  <big-figure|<image|../figures/module_2/mutation_genotype_prior.png|390pt|178pt||>|<label|fig:genotype_prior>Schematic
  of how to illicit mutational genotype priors. We assume we have the
  information from Figure <reference|fig:cn_snv> and know the major copy
  number is 2 and the minor copy number is 1. The first two examples
  correspond to mutations which happen prior to the copy number change, hence
  the total copy number of the reference and variant population differ. The
  third example corresponds to the case where the mutation occurs after the
  copy number event. Hence, the copy number of the reference and variant
  population are the same. Furthermore, only a single copy can be mutated by
  the infinite sites assumption.>

  Now, the question is how to use this copy number information to develop a
  prior for the mutational genotypes of the populations with respect to an
  SNV? There are many ways to do this, and it is really a matter of personal
  belief. The model we will adopt breaks into two cases.\ 

  <\enumerate-numeric>
    <item>If the mutation occurs before the copy number event, then we need
    to consider all possible mutational genotypes for the variant population
    with mutations on one to the major copy number of the chromosomes. Since
    the copy number event does not affect the reference population, the total
    copy number of the reference population will be the same as the normal
    population.

    <item>If the mutation occurs after the copy number event, then only a
    single chromosome can be mutated. This follows from the infinite sites
    assumption. We also have that the total copy number of the reference
    population matches the variant population.
  </enumerate-numeric>

  Applying these two rules we can list a set of possible mutational genotypes
  compatible with the observed copy number profile. We also need a way to
  assign probabilities to each of these scenarios. In the absence of any
  reason to prefer one over the other, we simply weight them equally.

  <subsubsection|The full independent model>

  Now we have a way to illicit priors for the mutational genotypes, we will
  incorporate this information into the model. Let <math|\<b-pi\>> be a
  vector of probabilities for each mutational genotype. We will use the
  notation <math|\<pi\><rsub|\<b-psi\>>> to indicate the prior probability of
  the mutational genotype <math|\<b-psi\>>. We can then compute the
  probability of the observed data as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|B=b,\<b-pi\>,\<phi\>,t,d|)>>|<cell|=>|<cell|p<around*|(|\<phi\>|)><big|sum><rsub|\<b-psi\>>p<around*|(|\<b-psi\>\|\<b-pi\>|)>p<around*|(|B=b\|\<b-psi\>,\<phi\>,t,d|)>>>|<row|<cell|>|<cell|=>|<cell|p<around*|(|\<phi\>|)><big|sum><rsub|\<b-psi\>>\<pi\><rsub|\<b-psi\>>
    p<around*|(|B=b\|\<b-psi\>,\<phi\>,t,d|)>>>>>
  </eqnarray>

  Here we have applied the law of total probability to sum over all possible
  genotypes we believe possible. Also note that
  <math|p<around*|(|B=b\|\<b-psi\>,\<phi\>,t,d|)>> is just the conditional
  distribution we defined earlier. Thus our full model is now

  <\eqnarray>
    <tformat|<table|<row|<cell|\<phi\>>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|0,1|)>>>|<row|<cell|B\|\<b-pi\>,\<phi\>,t,d>|<cell|\<sim\>>|<cell|<big|sum><rsub|\<b-psi\>>\<pi\><rsub|\<b-psi\>>
    <text|Binomial><around*|(|\<cdot\>\|d,\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)>>>>>
  </eqnarray>

  where the last line is a mixture distribution over the mutational
  genotypes.

  <big-figure|<image|../figures/module_2/density.pdf|424pt|279pt||>|<label|fig:ccf_density>Example
  posterior densities for the cellular prevalence <math|\<phi\>> computed
  from the model. The top row shows the case for a homozygous diploid
  position (there is a s typo it should be CN=(2,0)). The second row shows
  the posterior for the case illustrated in Figures <reference|fig:cn_snv>
  and <reference|fig:genotype_prior>.>

  We are now at a point where we can compute the posterior of <math|\<phi\>>
  by applying Bayes' rule. We cannot analytically compute the integral to get
  the normalisation constants, so we use numerical approximations. The result
  of doing this is illustrated in Figure <reference|fig:ccf_density>. One
  thing to note is that the posterior distribution is multi-modal. This
  corresponds to our uncertainty about the true genotype. For example in the
  first row of Figure <reference|fig:ccf_density> the VAF is 0.5 and the mode
  at 0.5 corresponds to the case where both copies are mutated, while the
  second mode at 1.0 corresponds to the case where a single copy is mutated.

  To recap, we have developed a model that allows us to infer the cellular
  prevalence of a mutation. We assume we observe allelic count data for an
  SNV (<math|a,b>), have an estimate of the tumour content (<math|t>) and
  know the copy number profile overlapping the SNV to derive the genotype
  prior (<math|\<b-pi\>>). This model treats all mutations as completely
  independent. A weakness of this approach is that the posteriors tend to be
  multi-modal, so we are still quite uncertain about the value of
  <math|\<phi\>>. In the next section we will discuss how to fix this
  problem.

  <subsection|Clustering mutations>

  In the previous section we systematically developed a model to infer the
  cellular prevalence of a single mutation. Using this model we can easily
  compute the posterior distribution for the cellular prevalence
  <math|\<phi\>> of any mutation. However, these posteriors tend to be
  multi-modal because we have to consider a large number of mutational
  genotypes. We will now look at how clustering mutations can help solve this
  problem.

  <subsubsection|Motivation>

  <big-figure|<image|../figures/module_2/phylogeny.pdf|261pt|105pt||>|<label|fig:snv_phylogeny>Illustration
  of the relationship between evolutionary history and cellular prevalence.
  On the left we have hypothetical evolutionary history, where stars indicate
  mutations and nodes clonal populations. On the write is a hypothetical set
  of cellular prevalence for the mutations.>

  Cancer is an evolutionary process and clonal populations are related by a
  phylogenetic tree. One important implication of this is that mutations
  which share the same evolutionary history will be at the same cellular
  prevalence. We illustrate this in Figure <reference|fig:snv_phylogeny>
  where the phylogeny is shown on the left and the cellular prevalence of
  mutations on the right. The most interesting mutations are the green, blue
  and purple which all appear at the some point in the phylogeny. Assuming
  mutations are propagated to descendants (never lost) then these mutations
  will always appear in the same set of cells. Hence, their cellular
  prevalence will be identical. This has two important implication:

  <\enumerate-numeric>
    <item>We should expect sets of mutations to have the same cellular.

    <item>If we can identify which mutations have the same cellular
    prevalence we can infer which ones share the same evolutionary history.
  </enumerate-numeric>

  The first point tells us that we should not treat mutations independently
  in our model. It also means that we can share statistical strength between
  mutations and potentially reduce the uncertainty in our posterior
  distributions for <math|\<phi\>>. The problem we face is that we do not
  know which mutations belong together or what their cellular prevalence are.
  We also do not know how many clones are in the sample. We will see how to
  address this in the remainder of this module.

  <\remark>
    Lost mutations are not actually a problem. If mutations originate at the
    same point and are lost at the same point they will form a cluster of
    there own. This will lead to the inference of an additional clone, but
    not impact the cellular prevalence estimates.\ 

    If we want to reconstruct the phylogeny based on the assumption mutations
    at higher prevalence are further up the tree, then we have a problem. We
    discuss this in module 4.
  </remark>

  <subsubsection|Mixture models>

  Before diving into the full model, we first review mixture models. Mixture
  models are a very useful type of probabilistic model which posit that
  datapoints originate from groups or clusters. Datapoints from the same
  cluster share the same parameters, and thus are similar to each other in
  some sense. For now we assume the number of clusters, <math|K>, is known in
  advance and fixed. Then a standard way to construct a Bayesian mixture
  model is as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-kappa\>>|<cell|\<in\>>|<cell|\<bbb-R\><rsub|+><rsup|K>>>|<row|<cell|\<b-rho\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|z<rsub|n>\|\<b-rho\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-rho\>|)>>>|<row|<cell|\<theta\><rsub|k>>|<cell|\<sim\>>|<cell|G<around*|(|\<cdot\>|)>>>|<row|<cell|x<rsub|n>\|\<b-theta\>,z<rsub|n>>|<cell|\<sim\>>|<cell|F<around*|(|\<cdot\>\|\<theta\><rsub|z<rsub|n>>|)>>>>>
  </eqnarray>

  This model associates a latent variable, <math|z<rsub|n>>, with each data
  point. The variable <math|z<rsub|n>> takes values in the set
  <math|<around*|{|1,\<ldots\>,K|}>>, and acts as an indicator for which
  cluster a data point originates from. Each cluster has an associated
  parameter <math|\<theta\><rsub|k>> sampled independently from a
  distribution <math|G>. The observed data <math|X<rsub|n>> is then generated
  from a distribution <math|F> with parameter
  <math|\<theta\><rsub|z<rsub|n>>>. Thus whenever
  <math|z<rsub|i>=z<rsub|j>=k> for data points <math|i> and <math|j> they are
  in the same cluster and have been generated from the same distribution
  <math|F<around*|(|\<theta\><rsub|k>|)>>.

  <subsubsection|Sharing statistical strength>

  We will now improve the previous model by turning it into a mixture model.
  For now we assume the number of clusters <math|K> is known and fixed. Then
  the updated model is\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-rho\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|z<rsub|n>\|\<b-rho\>>|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|\<b-rho\>|)>>>|<row|<cell|\<phi\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>>|<row|<cell|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>>|<cell|\<sim\>>|<cell|<big|sum><rsub|\<b-psi\>>\<pi\><rsub|n\<b-psi\>>
    <text|Binomial><around*|(|\<cdot\>\|d<rsub|n>,\<xi\><around*|(|\<b-psi\>,\<phi\><rsub|z<rsub|n>>,t|)>|)>>>>>
  </eqnarray>

  where we introduce the index <math|n> for data points. In this model data
  points are no longer independent, but will share the same cellular
  prevalence when they belong to the same cluster.

  This example nicely illustrates the modularity of Bayesian probabilistic
  models. Specifically, we were able to reuse the previous model for
  mutational genotypes and embed it in a more complex model. This is a useful
  strategy in general. Begin by breaking down the problem into simpler
  sub-problems, and then progressively extend the model.

  <subsubsection|MCMC inference>

  To fit the new model to the data, we can no longer appeal to numerical
  methods to compute the posterior because we have many more parameters, some
  of which are discrete. To address this we will use MCMC to approximate the
  posterior.\ 

  <\remark>
    We could quite easily use the expectation maximisation (EM) algorithm to
    compute the MAP estimate for this model rather than MCMC. For finite
    mixture models this will often be faster. The main disadvantage is that
    we will not have a full posterior, but point estimates as discussed in
    module 1. Thus we cannot quantify uncertainty. Once we move to using a
    Dirichlet process prior, EM will no longer be applicable so we will have
    to use MCMC.
  </remark>

  Before we work out the details we will write down the joint distribution.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-b\>,\<b-d\>,\<b-pi\>,\<b-phi\>,\<b-rho\>,\<b-kappa\>,t,\<b-z\>|)>>|<cell|=>|<cell|p<around*|(|\<b-rho\>\|\<b-kappa\>|)><big|prod><rsub|k=1><rsup|K>p<around*|(|\<phi\><rsub|k>|)><big|prod><rsub|n=1><rsup|N>p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>|)>
    p<around*|(|z<rsub|n>\|\<b-rho\>|)>>>>>
  </eqnarray>

  We will update the model parameters in blocks. This is typically done in
  MCMC methods, as designing good updates for all the parameters
  simultaneously is usually hard (see module 1). We will use a combination of
  Metropolis-Hastings (MH) and Gibbs sampling. The updates we will use are:

  <\itemize-dot>
    <item><math|\<b-rho\>> we will use a Gibbs update

    <item><math|\<phi\><rsub|k>> we will use an MH update with a
    <math|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>> proposal

    <item><math|z<rsub|n>> we will use a Gibbs update
  </itemize-dot>

  \;

  To implement the Gibbs step for <math|z<rsub|n>> we need to compute the
  conditional distribution <math|p<around*|(|z<rsub|n>\|-|)>>, where <math|->
  indicates all the other variables. We introduce the notation
  <math|\<b-z\><rsup|<around*|(|-n|)>>=<around*|(|z<rsub|1>,\<ldots\>,z<rsub|n-1>,z<rsub|n+1>,\<ldots\>,z<rsub|N>|)>>
  which is the vector of all cluster indicator variables except the
  <math|n<rsup|th>> one. The conditional distribution is then

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|z<rsub|n>=k\|-|)>>|<cell|=>|<cell|<frac|p<around*|(|\<b-b\>,\<b-d\>,\<b-pi\>,\<b-phi\>,\<b-rho\>,\<b-kappa\>,t,\<b-z\><rsup|<around*|(|-n|)>>,z<rsub|n>=k|)>|p<around*|(|\<b-b\>,\<b-d\>,\<b-pi\>,\<b-phi\>,\<b-rho\>,\<b-kappa\>,t,\<b-z\><rsup|<around*|(|-n|)>>|)>>>>|<row|<cell|>|<cell|=>|<cell|<frac|p<around*|(|\<b-b\>,\<b-d\>,\<b-pi\>,\<b-phi\>,\<b-rho\>,\<b-kappa\>,t,\<b-z\><rsup|<around*|(|-n|)>>,z<rsub|n>=k|)>|<big|sum><rsub|\<ell\>=1><rsup|K>p<around*|(|\<b-b\>,\<b-d\>,\<b-pi\>,\<b-phi\>,\<b-rho\>,\<b-kappa\>,t,\<b-z\><rsup|<around*|(|-n|)>>,z<rsub|n>=\<ell\>|)>>>>|<row|<cell|after
    some cancellation>|<cell|=>|<cell|<frac|\<rho\><rsub|k>
    p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>=k|)>|<big|sum><rsub|\<ell\>=1><rsup|K>\<rho\><rsub|\<ell\>>
    p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>|)>>>>>>
  </eqnarray>

  let <math|<wide|\<rho\>|\<bar\>><rsub|n k>=<frac|\<rho\><rsub|k>
  p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>=k|)>|<big|sum><rsub|\<ell\>=1><rsup|K>\<rho\><rsub|\<ell\>>
  p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<b-phi\>,t,d<rsub|n>,z<rsub|n>|)>>>
  then

  <\eqnarray>
    <tformat|<table|<row|<cell|z<rsub|n>\|->|<cell|\<sim\>>|<cell|<text|Categorical><around*|(|\<cdot\>\|<wide|\<b-rho\>|\<bar\>><rsub|n>|)>>>>>
  </eqnarray>

  \;

  Because we choose a Dirichlet prior for <math|\<b-rho\>> the conditional
  distribution is easily obtained from conjugacy, and will be a Dirichlet
  distribution as well. Let <math|m<rsub|k>=<big|sum><rsub|n=1><rsup|N>\<bbb-I\><around*|(|z<rsub|n>=k|)>>
  be the number of data points from cluster <math|k>. Let
  <math|<wide|\<b-kappa\>|\<bar\>>=<around*|(|\<kappa\><rsub|1>+m<rsub|1>,\<ldots\>,\<kappa\><rsub|k>+m<rsub|k>,\<ldots\>,\<kappa\><rsub|K>+m<rsub|K>|)>>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-rho\><mid|\|>->|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|<wide|\<b-kappa\>|\<bar\>>|)>>>>>
  </eqnarray>

  <subsubsection|Dirichlet process>

  The finite mixture model we defined previously has one major problem: We
  assume the number of clones <math|K> is known in advance. In practice this
  is not true, and we would like to infer the number of clones as part of the
  model. One strategy to address this problem is to use a Dirichlet process
  prior for the cellular prevalence.

  The Dirichlet process (DP) prior is an example of a
  <with|font-shape|italic|non-parametric> Bayesian prior. Informally these
  means that it is a prior which can adapt model complexity as more data is
  observed. More formally the DP is a distribution over distribution
  (stochastic process). This means that the random variables we sample from a
  DP are distributions. There are two parameters for the DP: the
  concentration parameter <math|\<alpha\>\<in\>\<bbb-R\><rsub|+>> \ and the
  base measure <math|G<rsub|0>> a distribution. Roughly speaking
  <math|\<alpha\>> controls how many clusters we expect. While the
  distribution <math|G<rsub|0>> is used to sample the new values that the
  distribution exhibits.

  To understand how a DP can be useful we need to take a slightly different
  view of mixture models. So far we have used the cluster indicator variables
  <math|z<rsub|n>> to identify which data points belong to the same cluster.
  So if <math|z<rsub|i>=z<rsub|j>=k> that means data points <math|i> and
  <math|j> come from a distribution with parameter <math|\<theta\><rsub|k>>,
  in other words belong to cluster <math|k>. We can change our viewpoint
  though, and instead of having one <math|\<theta\><rsub|k>> for each
  cluster, we can instead assign data point <math|n> its own parameter
  <math|\<zeta\><rsub|n>>. When we have <math|\<zeta\><rsub|i>=\<zeta\><rsub|j>>
  we then say data points <math|i> and <math|j> belong to the sample cluster.
  Now if we sample <math|\<zeta\><rsub|n>> from a continuous distribution,
  then there is zero probability that two data points will ever share the
  same value. So we need to sample <math|\<zeta\><rsub|n>> from a discrete
  distribution if there is any chance for shared values. In the background
  this is what the mixture model is doing, it is creating the discrete
  distribution

  <\eqnarray>
    <tformat|<table|<row|<cell|G>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|K>\<rho\><rsub|k>
    \<delta\><rsub|\<theta\><rsub|k>><around*|(|\<cdot\>|)>>>>>
  </eqnarray>

  where <math|\<delta\><rsub|x><around*|(|y|)>> is the Kronecker delta
  function that equals one when <math|x=y> or zero otherwise. We then draw
  <math|\<zeta\><rsub|n>\<sim\>G<around*|(|\<cdot\>|)>>, so that
  <math|\<zeta\><rsub|n>\<in\><around*|{|\<theta\><rsub|1>,\<ldots\>,\<theta\><rsub|K>|}>>.\ 

  We can re-write our current mixture model with this viewpoint as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-rho\>\|\<b-kappa\>>|<cell|\<sim\>>|<cell|<text|Dirichlet><around*|(|\<cdot\>\|\<b-kappa\>|)>>>|<row|<cell|\<theta\><rsub|k>>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>>|<row|<cell|G>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|K>\<rho\><rsub|k>
    \<delta\><rsub|\<theta\><rsub|k>><around*|(|\<cdot\>|)>>>|<row|<cell|\<phi\><rsub|n>>|<cell|\<sim\>>|<cell|G>>|<row|<cell|b<rsub|n>\|\<b-pi\><rsub|n>,\<phi\><rsub|n>,t,d<rsub|n>,z<rsub|n>>|<cell|\<sim\>>|<cell|<big|sum><rsub|\<b-psi\>>\<pi\><rsub|n\<b-psi\>>
    <text|Binomial><around*|(|\<cdot\>\|d<rsub|n>,\<xi\><around*|(|\<b-psi\>,\<phi\><rsub|n>,t|)>|)>>>>>
  </eqnarray>

  Note that the cluster indicators are no longer present and each data point
  no samples its own cellular prevalence <math|\<phi\><rsub|n>>. However, the
  values must belong to the set <math|<around*|{|\<theta\><rsub|1>,\<ldots\>,\<theta\><rsub|K>|}>>
  so we have clustering. Here the <math|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>
  is used to identify the set of values <math|G> takes on. The elements of
  this set are often referred to as the atoms of the distribution.

  The previous discussion may seem like a very confusing way to define a
  mixture model. The reason that it is useful is that if sample a
  distribution <math|G> from a DP it will be discrete. More explicitly, any
  distribution <math|G> sample from a DP has the following form

  <\eqnarray>
    <tformat|<table|<row|<cell|G>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|\<infty\>>\<rho\><rsub|k>
    \<delta\><rsub|\<theta\><rsub|k>><around*|(|\<cdot\>|)>>>>>
  </eqnarray>

  Thus we can use this distribution <math|G> in our mixture model, just the
  \ same way as the finite case. The updated version of our model now takes
  the following form.

  \;

  <\eqnarray>
    <tformat|<table|<row|<cell|G<rsub|0>>|<cell|=>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]>|)>>>|<row|<cell|G\|\<alpha\>,G<rsub|0>>|<cell|\<sim\>>|<cell|<text|DP><around*|(|\<cdot\>\|\<alpha\>,G<rsub|0>|)>>>|<row|<cell|\<phi\><rsub|n>>|<cell|\<sim\>>|<cell|G>>|<row|<cell|b<rsub|n>\|\<b-pi\><rsub|n>,\<phi\><rsub|n>,t,d<rsub|n>,z<rsub|n>>|<cell|\<sim\>>|<cell|<big|sum><rsub|\<b-psi\>>\<pi\><rsub|n\<b-psi\>>
    <text|Binomial><around*|(|\<cdot\>\|d<rsub|n>,\<xi\><around*|(|\<b-psi\>,\<phi\><rsub|n>,t|)>|)>>>>>
  </eqnarray>

  <subsubsection|Chinese restaurant process>

  Before discussin how to fit the model we will take a brief digression to
  discuss the Chinese restaurant process (CRP). The CRP is the marginal
  distribution of the DP when we integrate out the distribution <math|G>.
  Formally, the CRP is a probability distribution on partitions of the
  integers. Let <math|<around*|[|n|]>=<around*|{|1,\<ldots\>,n|}>> be the set
  of all positive integers up to <math|n>. Then a partition of <math|n> is a
  set <math|c<rsub|n>=<around*|{|b:b\<subset\><around*|[|n|]>|}>> such that
  <math|b\<cap\>b<rprime|'>=\<varnothing\>> and
  <math|<big|cup><rsub|b\<in\>c<rsub|n>>b=<around*|[|n|]>>. That is a set of
  disjoint sets whose union equals <math|<around*|[|n|]>>. With some thought
  you can see this is equivalent to a clustering of the data points labelled
  from <math|><math|1> to <math|n>. The probability mass function (pmf) of
  the CRP is given by

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|c<rsub|N>\|\<alpha\>|)>>|<cell|=>|<cell|<frac|\<Gamma\><around*|(|\<alpha\>|)>|\<Gamma\><around*|(|\<alpha\>+N|)>>\<alpha\><rsup|<around*|\||c|\|>>
    <big|prod><rsub|b\<in\>c<rsub|N>><around*|(|<around*|\||b|\|>-1|)>!>>>>
  </eqnarray>

  where <math|\<alpha\>> is the concentration parameter.

  The CRP is often described by an analogy to customers entering a chinese
  restaurant, hence the name. The description goes as follows. The first
  customer enters the restauran and sits down at a table. The second customer
  then enters the restaurant makes a choice. The can either join the first
  customer with probability <math|<frac|1|1+\<alpha\>>> or start a new table
  with probability <math|<frac|\<alpha\>|1+\<alpha\>>>. As new customers
  enter they can choose to set at an existing table with probability
  proportional to the number of customers already there, or they can start a
  new table with probability proportional to <math|\<alpha\>>. The
  distribution over seatings of customers is then given by the CRP pmf.

  There are a few interesting aspects of this process. First, it has a rich
  get richer property where new customers are more likely to join existing
  tables. Second, the process is <with|font-shape|italic|exchangeable> so the
  order customers enter the restaurant does not affect the distribution. This
  can be seen directly from the pmf of the CRP. This feature is particularly
  useful, as it will allow to develop a Gibbs sampler for the DP in the next
  section.

  <subsubsection|Full model inference>\ 

  We can still use MH updates for <math|\<phi\><rsub|n>>. The new problem is
  to update <math|\<phi\><rsub|n>>. There are two approaches to performing
  inferences. One is to make use of the <with|font-shape|italic|stick
  breaking representation> of the DP. Using this approach we can sample the
  distribution <math|G> directly. Once we have this distribution, then we can
  perform Gibbs updates just like the finite mixture model case. The second
  approach which we will follow is marginalise <math|G> and in which case we
  have a CRP.\ 

  The key insight is that we can use exchangeability to treat the data point
  we want to update as the last data point added. We will re-introduce the
  cluster labels for this step, let <math|z<rsub|n>\<in\><around*|{|1,\<ldots\>,K|}>>
  be the cluster label for the <math|n<rsup|th>> data point. Here <math|K> is
  the number of clusters with data points, after removing the
  <math|n<rsup|th>> data point. The the probability of this data point
  joining an existing cluster depends on the size of the cluster, excluding
  the data point of interest. We will use
  <math|m<rsub|k><rsup|<around*|(|-n|)>>=<big|sum><rsub|n>\<bbb-I\><around*|(|z<rsub|n>=k|)>>
  to denote these cluster sizes. The conditional probability of
  <math|z<rsub|n>> given the other indicator variables,
  <math|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>> is given by the following
  equation.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|z<rsub|n>=k\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|<frac|m<rsub|k><rsup|<around*|(|-n|)>>|n-1+\<alpha\>>>|<cell|<text|if>>|<cell|k\<in\><around*|{|1,\<ldots\>,K|}>>>|<row|<cell|<frac|\<alpha\>|n-1+\<alpha\>>>|<cell|<text|if>>|<cell|k=K+1>>>>>>>>>
  </eqnarray>

  <\remark>
    This update is fairly easy to understand. The only subtlety when we are
    updating a data point that is in a singleton cluster (only member of the
    cluster). In this case <math|K> differs when remove the data point from
    the clustering.
  </remark>

  This provides the conditional probability for the cluster label, but we
  will actually need a conditional probability of the form
  <math|p<around*|(|z<rsub|n>\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>,<around*|{|\<theta\><rsub|k>|}><rsub|k=1><rsup|K>,X|)>>
  since there is a cluster parameter associdated with each cluster and this
  is used to generate the data point.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|z<rsub|n>=k\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>,\<theta\>,X|)>>|<cell|=>|<cell|<frac|p<around*|(|X\|z<rsub|n>=k,<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>,\<theta\>|)>p<around*|(|z<rsub|n>=k\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>|)>|p<around*|(|X\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>,\<theta\>|)>>>>|<row|<cell|>|<cell|\<propto\>>|<cell|p<around*|(|x<rsub|n>\|z<rsub|n>=k,\<theta\>|)>
    p<around*|(|z<rsub|n>=k\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>|)>>>>>
  </eqnarray>

  Here <math|p<around*|(|x<rsub|n>\|z<rsub|n>=k,\<theta\>|)>> would just be
  the data likelihood for data point <math|n> when it has parameter
  <math|\<theta\><rsub|k>>. One tricky issue is that we do not have
  parameters for the new clusters, so we cannot easily evaluate this. One
  common approach is to integrate out the model parameters. Unfortunately
  this only works if the likelihood is conjugate to the prior, which is not
  the case for us. There are few ways around this. The one we will employ is
  to draw <math|l> values of <math|\<theta\>> from the prior to create
  <math|l> possible empty clusters to join. There is one last subtle point.
  When <math|n> was originally part of the singleton cluster we need to use
  that value as one of the <math|l> values of <math|\<theta\>>, so we only
  sample <math|l-1> new values from the prior. The probability of joining
  these clusters is then modified and becomes
  <math|<frac|<frac|\<alpha\>|l>|n-1+\<alpha\>>>. The number of empty table
  paramters <math|l> is a tuning paramter for the algorithm. In practice we
  typically take <math|l=2>. The Gibbs update is then given by the following
  formula.\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|z<rsub|n>=k\|<around*|{|z<rsub|i>|}><rsub|i\<neq\>n>,\<theta\>,X|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|<frac|m<rsub|k><rsup|<around*|(|-n|)>>|n-1+\<alpha\>>
    p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<phi\><rsub|k>,t,d<rsub|n>|)>
    >|<cell|<text|if>>|<cell|k\<in\><around*|{|1,\<ldots\>,K|}>>>|<row|<cell|<frac|<frac|\<alpha\>|l>|n-1+\<alpha\>>
    p<around*|(|b<rsub|n>\|\<b-pi\><rsub|n>,\<phi\><rsub|k>,t,d<rsub|n>|)>>|<cell|<text|if>>|<cell|k=K+1>>>>>>>>>
  </eqnarray>

  <\remark>
    The trick for sampling parameters for the new clusters may seem a bit
    obscure. For more details see [cite Neal] for this (algorithm 8 from the
    paper) and other possible ways to peform inference for DPs.\ 
  </remark>

  It is also useful to use an MH step to update the values
  <math|\<phi\><rsub|k>> between the updates of <math|z<rsub|n>>. This will
  allow us to move from values sampled from the prior, to those that are a
  close fit to the data. Without this MH move we would only update the
  cellular prevalences when we sample new tables. This is extremely slow we
  do not use any information about the data to generate the new
  <math|\<phi\>> values.

  <subsubsection|Computing consensus clustering>

  Thus far we have defined the model and come up with a strategy to fit the
  model. This involves running an MCMC algorithm for many iterations and
  collecting samples. At each iteration we will record the cluster paramters
  (cellular prevalences) and the cluster indicators (which cluster a
  datapoint was assigned to). Each iteration of the MCMC will have a
  different clustering of the data. The question is then how pick a best
  clustering?\ 

  One simple approach would be to take the sample with the highest joint
  probability i.e. the MAP estimator. This is sub-optimal as we ignore all
  the other samples from the MCMC chain. A better approach is to use a loss
  function. We are free to choose any loss function, but there are some
  challenges we need to consider. The biggest issue is that we can permute
  the labels of the cluster indicators and the cluster parameters and the
  likelihood is the same. So though a data point may move from cluster
  <math|k> to <math|\<ell\>>, nothing has changed because we also changing
  the cluster paramters from <math|\<theta\><rsub|k>> to
  <math|\<theta\><rsub|\<ell\>>>. Thus we ideally want a loss function that
  is invariant to label permutations. The other issue is that the number of
  clusters varies between iterations.

  Here we will use the adjusted rand index (ARI) as the loss function. This
  is a measure of clustering similarity. We will then seek the clustering
  which minimises this lost under our approximate posterior. The details on
  how to do this can be found in [cite mpear]. To implement this proceedure
  we compute the pair-wise similarity matrix of two data points. This is the
  proportion of MCMC samples in which that data points belong to the same
  cluster. A nice feature of this summary of the MCMC trace is that it does
  not depend on the number of clusters or the actual labels of the data
  points. To optimise the ARI we then build a dendrogram from the similarity
  matrix and find the cut level which maximises the MPEAR score defined in
  [cite mpear]. This yields our estimated clustering of the data
  <math|<wide|\<b-z\>|^>=<around*|(|<wide|z|^><rsub|1>,\<ldots\>,<wide|z|^><rsub|N>|)>>.

  We would also like to obtain posterior distributions for the cellular
  prevalence. Once approach is to look at the full MCMC trace of the cellular
  prevalence associated with a mutation. We can plot these values as
  histogram or using a density estimator. We can also report the mean and
  variance of this distribution. One downside of this approach is we report a
  separate posterior for each mutation, which ignores our estimate clustering
  <math|<wide|\<b-z\>|^>>. If we want to get the posterior cellular
  prevalence give the clustering we can use a slightly ad-hoc approach and
  compute the conditional posterior given <math|<wide|\<b-z\>|^>>. Suppose we
  want the conditional posterior for cluster <math|k>, that is the posterior
  of <math|\<phi\><rsub|k>>. Then we need to computer the following quantity.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<phi\><rsub|k>\|<wide|\<b-z\>|^>,X|)>>|<cell|=>|<cell|<frac|p<around*|(|\<phi\><rsub|k>|)><big|prod><rsub|<around*|{|n:z<rsub|n>=k|}>>p<around*|(|x<rsub|n>\|\<phi\><rsub|k>|)>|<big|int>p<around*|(|\<phi\><rsub|k>|)><big|prod><rsub|<around*|{|n:z<rsub|n>=k|}>>p<around*|(|x<rsub|n>\|\<phi\><rsub|k>|)>
    \<mathd\>\<phi\><rsub|k>>>>>>
  </eqnarray>

  The integral in the denomintor does not have a closed form. But it is a one
  dimensional integral so we can use numerical methods to compute it.

  <subsubsection|Multiple samples>

  Thus far we have only considered data from a single sample. If multiple
  samples are available it would be useful to model them jointly so that
  clusters are coherently defined across samples. Multiple samples are very
  useful for this problem, as they let us identify clusters which may be
  missed from single sample analysis. One reason this could happen is that
  some clones are only present in a subset of samples. Another reason is that
  we may have two clones at similar prevalence in one sample that would be
  impossible to identify. If we have another sample where there prevalence is
  quite different, then we have a better chance of finding both clones.

  The changes required to include multiple samples are actually fairly
  simple. We will need some notation for this. Let <math|m> be the index over
  samples, <math|M> be the number of samples, <math|b<rsub|n m>> indicate the
  number of reads with variant <math|n> in sample <math|m>, <math|d<rsub|n
  m>> be the corresponding read depth, <math|\<b-pi\><rsub|n m>> be the
  genotype prior for mutation <math|n> in sample <math|m> and
  <math|t<rsub|m>> be the tumour content of the <math|m<rsup|th>> sample.
  Then the new model becomes

  \;

  <\eqnarray>
    <tformat|<table|<row|<cell|G<rsub|0>>|<cell|=>|<cell|<text|Uniform><around*|(|\<cdot\>\|<around*|[|0,1|]><rsup|M>|)>>>|<row|<cell|G\|\<alpha\>,G<rsub|0>>|<cell|\<sim\>>|<cell|<text|DP><around*|(|\<cdot\>\|\<alpha\>,G<rsub|0>|)>>>|<row|<cell|\<b-phi\><rsub|n>>|<cell|\<sim\>>|<cell|G>>|<row|<cell|b<rsub|n
    m>\|\<b-pi\><rsub|n m>,\<phi\><rsub|n>,t<rsub|m>,d<rsub|n
    m>,z<rsub|n>>|<cell|\<sim\>>|<cell|<big|sum><rsub|\<b-psi\>>\<pi\><rsub|n
    m\<b-psi\>> <text|Binomial><around*|(|\<cdot\>\|d<rsub|n
    m>,\<xi\><around*|(|\<b-psi\>,\<phi\><rsub|n m>,t<rsub|m>|)>|)>>>>>
  </eqnarray>

  The main change is that we use Uniform prior over the <math|M> dimensional
  unit cube and we sample a vector of cellular prevalence
  <math|\<b-phi\><rsub|n>=<around*|(|\<phi\><rsub|n
  1>,\<ldots\>,\<phi\><rsub|n m>|)>> for each mutation. The inference
  procedure is largely the same. The only major differences is that we will
  now need to compute summaries of cellular prevalence for each sample.

  We see again the modularity of the Bayesian modelling approach. It is a
  fairly trivial exercise to extend our simpler model to the more complex
  case. The most tedious part is updating to a decent notation.

  <subsubsection|Overdispersion>

  One issue with current model is that read counts often except more
  variability than the Binomial can model. The issue is referred to as
  overdispersion. At the depths commonly used for WGS (30x-100x) this problem
  is usually not obvious. In higher coverage data such as targeted sequencing
  (<math|10<rsup|3>>x-<math|10<rsup|5>>x) it becomes more pronounced.\ 

  The solution is to use a distribution with more parameters and more
  flexibility. In the case of the Binomial, a common choice is to use the
  Beta-Binomial distribution which is overdispersed relative to the binomial.
  The Beta-Binomial has to parameters <math|a> and <math|b> like a Beta
  distribution. We can reparameterise the Beta-Binomial in terms of mean and
  variance as well. Using this approach we can set the mean to <math|\<xi\>>
  and fit the variance parameter. The procedure of substituting a more
  flexible distribution with the same mean is a common approach.

  <subsection|Discussion>

  In this module we saw how we to construct a Bayesian probabilistic model
  for inferring clonal population structure from bulk data. We started with a
  simple model to correct for mutational genotype and normal contamination.
  We then considered a more complex mixture model, ultimately using a
  Dirichlet process to infer the number of clones. We discussed how to fit
  the model using MCMC and summarise the resulting posterior approximation.

  This module illustrates the basic technique and steps needed to construct a
  probabilsitic model and fit it. The most important concepts are:

  <\enumerate-numeric>
    <item>Clearly define the problem

    <item>Start by solving simpler sub-problems

    <item>Iteratively extend the model

    <item>Identify a suitable method for fitting the model

    <item>If using MCMC, identify a way to report summaries of the posterior
  </enumerate-numeric>
</body>

<\initial>
  <\collection>
    <associate|preamble|false>
    <associate|project-flag|true>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-10|<tuple|5|?>>
    <associate|auto-11|<tuple|1.2.3|?>>
    <associate|auto-12|<tuple|1.2.4|?>>
    <associate|auto-13|<tuple|6|?>>
    <associate|auto-14|<tuple|7|?>>
    <associate|auto-15|<tuple|1.2.5|?>>
    <associate|auto-16|<tuple|8|?>>
    <associate|auto-17|<tuple|1.3|?>>
    <associate|auto-18|<tuple|1.3.1|?>>
    <associate|auto-19|<tuple|9|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-20|<tuple|1.3.2|?>>
    <associate|auto-21|<tuple|1.3.3|?>>
    <associate|auto-22|<tuple|1.3.4|?>>
    <associate|auto-23|<tuple|1.3.5|?>>
    <associate|auto-24|<tuple|1.3.6|?>>
    <associate|auto-25|<tuple|1.3.7|?>>
    <associate|auto-26|<tuple|1.3.8|?>>
    <associate|auto-27|<tuple|1.3.9|?>>
    <associate|auto-28|<tuple|1.3.10|?>>
    <associate|auto-29|<tuple|1.4|?>>
    <associate|auto-3|<tuple|1|?>>
    <associate|auto-4|<tuple|1.2|?>>
    <associate|auto-5|<tuple|1.2.1|?>>
    <associate|auto-6|<tuple|2|?>>
    <associate|auto-7|<tuple|3|?>>
    <associate|auto-8|<tuple|1.2.2|?>>
    <associate|auto-9|<tuple|4|?>>
    <associate|fig:bulk_sequencing|<tuple|1|?>>
    <associate|fig:ccf_density|<tuple|8|?>>
    <associate|fig:ccf_vs_vaf|<tuple|2|?>>
    <associate|fig:cn_snv|<tuple|6|?>>
    <associate|fig:genotype_prior|<tuple|7|?>>
    <associate|fig:mut_genotype|<tuple|3|?>>
    <associate|fig:pop_structure|<tuple|4|?>>
    <associate|fig:pop_structure_simple|<tuple|5|?>>
    <associate|fig:snv_phylogeny|<tuple|9|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<surround|<hidden-binding|<tuple>|1>||Schematic of bulk
      sequencing for a tumour experiment. At the top we have the input cell
      population, where stars indicate mutations. In the middle we show the
      aligned reads obtained from performing bulk sequencing. Positions in
      reads are colour coded to match the mutations at the top. Note the
      proportion of reads with a variant is roughly similar to the proportion
      of input cells with the mutation. At the bottom we show the summarised
      read counts which we will use for modelling.>|<pageref|auto-3>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|2>||Example of a
      heterozygous diploid mutation showing why the variant allele frequency
      (VAF) is not the same proportion of cells harbouring the mutation
      (cellular prevalence). On the left we have the input population of
      cells which all have the mutation. On the right we have the observed
      sequence data where only half the reads (on average) have the
      mutation.>|<pageref|auto-6>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|3>||Effect of
      mutational genotype on observed VAF.>|<pageref|auto-7>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|4>||Illustration of the
      assumed population structure. Here all populations are defined with
      respect to a single mutation. The circular cells are non-malignant and
      the irregularly shaped ones are malignant.>|<pageref|auto-9>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|5>||Illustration of the
      simplified population structure. In contrast to Figure
      <reference|fig:pop_structure_simple> the mutational genotypes of all
      cells within the reference and variant populations are the same. Note
      the genotypes are different between the populations, as they must be
      since one population has the mutation and the other does
      not.>|<pageref|auto-10>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|6>||Illustration of
      allele specific copy number profile. The red line is the major copy
      number and the blue line is the minor copy number. The star indicates
      the location of an SNV. The major copy number of this SNV is 2 and the
      minor is 1.>|<pageref|auto-13>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|7>||Schematic of how to
      illicit mutational genotype priors. We assume we have the information
      from Figure <reference|fig:cn_snv> and know the major copy number is 2
      and the minor copy number is 1. The first two examples correspond to
      mutations which happen prior to the copy number change, hence the total
      copy number of the reference and variant population differ. The third
      example corresponds to the case where the mutation occurs after the
      copy number event. Hence, the copy number of the reference and variant
      population are the same. Furthermore, only a single copy can be mutated
      by the infinite sites assumption.>|<pageref|auto-14>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|8>||Example posterior
      densities for the cellular prevalence <with|mode|<quote|math>|\<phi\>>
      computed from the model. The top row shows the case for a homozygous
      diploid position (there is a s typo it should be CN=(2,0)). The second
      row shows the posterior for the case illustrated in Figures
      <reference|fig:cn_snv> and <reference|fig:genotype_prior>.>|<pageref|auto-16>>

      <tuple|normal|<surround|<hidden-binding|<tuple>|9>||Illustration of the
      relationship between evolutionary history and cellular prevalence. On
      the left we have hypothetical evolutionary history, where stars
      indicate mutations and nodes clonal populations. On the write is a
      hypothetical set of cellular prevalence for the
      mutations.>|<pageref|auto-19>>
    </associate>
    <\associate|toc>
      1.<space|2spc>Inferring clonal population structure from SNV data
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2.<space|2spc>Modelling mutation genotype
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.2.1.<space|2spc>Background
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.2.2.<space|2spc>Population structure
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|1.2.3.<space|2spc>Modelling counts
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|1.2.4.<space|2spc>Genotype priors
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|1.2.5.<space|2spc>The full independent
      model <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|1.3.<space|2spc>Clustering mutations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|2tab>|1.3.1.<space|2spc>Motivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|2tab>|1.3.2.<space|2spc>Mixture models
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|2tab>|1.3.3.<space|2spc>Sharing statistical
      strength <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|2tab>|1.3.4.<space|2spc>MCMC inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|2tab>|1.3.5.<space|2spc>Dirichlet process
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|2tab>|1.3.6.<space|2spc>Chinese restaurant
      process <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|2tab>|1.3.7.<space|2spc>Full model inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25>>

      <with|par-left|<quote|2tab>|1.3.8.<space|2spc>Computing consensus
      clustering <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <with|par-left|<quote|2tab>|1.3.9.<space|2spc>Multiple samples
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27>>

      <with|par-left|<quote|2tab>|1.3.10.<space|2spc>Overdispersion
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28>>

      <with|par-left|<quote|1tab>|1.4.<space|2spc>Discussion
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-29>>
    </associate>
  </collection>
</auxiliary>