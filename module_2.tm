<TeXmacs|1.99.8>

<style|<tuple|tmbook|british>>

<\body>
  <section|Inferring clonal population structure from SNV data>

  In this module we consider the problem of inferring clonal population
  structure using bulk sequencing data. We specifically consider the case of
  using SNVs a our clonal markers, differing the use of copy number and
  rearrangement breakpoints to later modules.

  <subsection|Problem>

  We consider the problems of how to use high throughput sequencing to infer
  the clonal population structure of a tumour. This problem is somewhat old
  now, but still remain relevant when looking at large cancer cohort datasets
  that are continuing to be generated. The basic ideas of deconvolution are
  also appearing in other areas so this module should provide some useful
  insight.

  Precisely we will use read count data from SNVs to

  <\enumerate-numeric>
    <item>Infer what proportion of cancer cells harbour a mutation

    <item>Infer what mutations share the same evolutionary history i.e.
    originate at the same time and are lost in the same subsets of clones
  </enumerate-numeric>

  <big-figure|<image|figures/bulk_sequencing_v2.pdf||400pt||>|<label|fig:bulk_sequencing>Schematic
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
  fact that not all mutation will heterozygous diploid events due to
  coincident copy number variation. Normal contamination refers to the fact
  we also sequence non-malignant cells in real tumours. The second part will
  be a mechanism to cluster the mutations. Here we will use the formalism of
  mixture models. There is one challenging issue, which is that we do not
  know how many clones there are in the sample(s). We will address this
  problem using the Dirichlet process.

  <subsection|Modelling mutation genotype>

  The first issue we tackle is the problem of mutational genotype. Our goal
  in this section is to define a probabilistic model that links the observed
  read count data with cellular prevalence. We can then apply the standard
  Bayesian machinery to compute estimates of cellular prevalence.

  <subsubsection|Background>

  We will make the assumption that a point mutation only occurs once at a
  locus during the evolutionary history of tumour. This is often referred to
  as the <with|font-shape|italic|infinite sites assumption>. This is
  assumption is motivated by the fact that the mutation rate is usually
  relatively low compared to the total size of the genome. It can break down
  at mutational hot spots, and this has been observed, for example when
  multiple substitutions are observed at a single locus. So this assumption
  represents our first approximation to truth.

  <big-figure|<image|figures/allelic_vs_cellular_frequency.v2.pdf|461pt|219pt||>|<label|fig:ccf_vs_vaf>Example
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

  <big-figure|<image|figures/genotype.pdf|434pt|239pt||>|<label|fig:mut_genotype>Effect
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
    truly sampled the reference and observe an erroneous variant, or vice
    versa. Because of this symmetry the error term only appears in cases when
    all of the chromosomes have the reference or the variant allele.
  </enumerate-numeric>

  <big-figure|<image|figures/population_structure.eps|400pt|||>|<label|fig:pop_structure>Illustration
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

  <big-figure|<image|figures/population_structure_simple.eps|400pt|||>|<label|fig:pop_structure_simple>Illustration
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
    the mutation for genotype <math|g>. Me actually make a small modification
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
    If maybe surprising that we consider how many copies of the locus a cell
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
    which chromosomes from the cell get sequenced.
  </remark>

  Note we need to normalise the probabilities to sum to one. Formally, let
  <math|E<rsub|i>\<in\><around*|{|N,R,V|}>> be a random variable indicating
  that read <math|i> was sample from population <math|a>. Then

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|E<rsub|i>=e|)>>|<cell|=>|<cell|<choice|<tformat|<table|<row|<cell|<frac|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>|Z>>|<cell|<text|if>>|<cell|e=N>>|<row|<cell|<frac|t
    <around*|(|1-\<phi\>|)> c<around*|(|g<rsub|R>|)>|Z>>|<cell|<text|if>>|<cell|e=R>>|<row|<cell|<frac|t
    \<phi\> c<around*|(|g<rsub|V>|)>|Z>>|<cell|<text|if>>|<cell|e=V>>>>>>>|<row|<cell|Z>|<cell|=>|<cell|<around*|(|1-t|)>c<around*|(|g<rsub|N>|)>+t
    <around*|(|1-\<phi\>|)> c<around*|(|g<rsub|R>|)>+t \<phi\>
    c<around*|(|g<rsub|V>|)>>>>>
  </eqnarray>

  Once we know which population a read comes from, then it is straightforward
  to compute the probability the read has the mutation. This is simply given
  by <math|\<mu\><around*|(|g<rsub|e>|)>> if <math|E=e>. To formalise this
  let <math|F<rsub|i>> be a random variable indicating if read <math|i> is
  has the mutation. With this notation we have that

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
    <tformat|<table|<row|<cell|\<phi\>>|<cell|\<sim\>>|<cell|<text|Uniform><around*|(|\<cdot\>\|0,1|)>>>|<row|<cell|B\|\<b-psi\>,\<phi\>,t,d>|<cell|\<sim\>>|<cell|<text|Binomial><around*|(|\<cdot\>\|d,\<xi\><around*|(|\<b-psi\>,\<phi\>,t|)>|)>>>>>
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

  \;

  \;

  \;

  \;

  \;

  \;

  \;

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
    <associate|auto-1|<tuple|1|1|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-10|<tuple|5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-11|<tuple|1.2.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-2|<tuple|1.1|1|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-3|<tuple|1|1|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-4|<tuple|1.2|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-5|<tuple|1.2.1|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-6|<tuple|2|3|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-7|<tuple|3|4|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-8|<tuple|1.2.2|5|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-9|<tuple|4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:bulk_sequencing|<tuple|1|1|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:ccf_vs_vaf|<tuple|2|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:mut_genotype|<tuple|3|3|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:pop_structure|<tuple|4|4|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|fig:pop_structure_simple|<tuple|5|5|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<surround|<hidden|<tuple>>||Schematic of bulk sequencing
      for a tumour experiment. At the top we have the input cell population,
      where stars indicate mutations. In the middle we show the aligned reads
      obtained from performing bulk sequencing. Positions in reads are colour
      coded to match the mutations at the top. Note the proportion of reads
      with a variant is roughly similar to the proportion of input cells with
      the mutation. At the bottom we show the summarised read counts which we
      will use for modelling.>|<pageref|auto-3>>

      <tuple|normal|<surround|<hidden|<tuple>>||Example of a heterozygous
      diploid mutation showing why the variant allele frequency (VAF) is not
      the same proportion of cells harbouring the mutation (cellular
      prevalence). On the left we have the input population of cells which
      all have the mutation. On the right we have the observed sequence data
      where only half the reads (on average) have the
      mutation.>|<pageref|auto-5>>

      <tuple|normal|<surround|<hidden|<tuple>>||Effect of mutational genotype
      on observed VAF.>|<pageref|auto-6>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of the assumed
      population structure. Here all populations are defined with respect to
      a single mutation. The circular cells are non-malignant and the
      irregularly shaped ones are malignant.>|<pageref|auto-7>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of the
      simplified population structure. In contrast to Figure
      <reference|fig:pop_structure_simple> the mutational genotypes of all
      cells within the reference and variant populations are the same. Note
      the genotypes are different between the populations, as they must be
      since one population has the mutation and the other does
      not.>|<pageref|auto-8>>
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
    </associate>
  </collection>
</auxiliary>