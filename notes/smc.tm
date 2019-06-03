<TeXmacs|1.99.8>

<style|<tuple|tmbook|british>>

<\body>
  <section|Miscallaneous material>

  <subsection|Sequential Monte Carlo>

  Sequential Monte Carlo (SMC) is a powerful MCMC method which can be useful
  for sampling high dimensional variables. Classically SMC has been used for
  models with an obvious sequential structure such as HMMs. However, the
  method is more general and can be applied any time we can induce a
  sequential structure. In this section we will review the details of SMC.
  Before we introducing the full SMC algorithm we will look at some simpler
  sampling approaches which were precursors to SMC.

  <subsubsection|Important sampling>

  The basic problem we have is that we would like to compute Monte Carlo
  estimates to approximate the expectation of some function <math|g> with
  respect to some probability distribution density <math|p>. We will abuse
  notation and also use the density function when referring to the
  distribution. Assuming we can sample from <math|p> this means we would like
  to draw <math|S> values <math|<around*|{|x<rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>>
  to make the following approximation.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|p><around*|[|h|]>>|<cell|\<assign\>>|<cell|<big|int>h<around*|(|x|)>p<around*|(|x|)>\<mathd\>x>>|<row|<cell|>|<cell|\<approx\>>|<cell|<frac|1|S><big|sum><rsub|s=1><rsup|S>h<around*|(|x<rsup|<around*|(|s|)>>|)>>>>>
  </eqnarray>

  The connection with Bayesian inference is that <math|p<around*|(|x|)>> will
  be our posterior and <math|x> the model paramters. In practice sampling
  from <math|p> is often hard which is why we revert to MCMC methods.

  Importance sampling (IS) can be useful when it is possible to evaluate
  <math|p> but not sample from it. The basic idea is to introduce another
  distribution <math|q> which is easy to sample from and it possible to
  evaluate its density. When can then make use of the following simple
  identity.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|p><around*|[|h|]>>|<cell|=>|<cell|<big|int>h<around*|(|x|)><wide*|<frac|p<around*|(|x|)>|q<around*|(|x|)>>|\<wide-underbrace\>><rsub|w<around*|(|x|)>>
    q<around*|(|x|)>\<mathd\>x>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q><around*|[|h<around*|(|x|)>
    w<around*|(|x|)>|]>>>>>
  </eqnarray>

  Then we have the following Monte Carlo approximation

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|p><around*|[|h|]>>|<cell|\<approx\>>|<cell|<frac|1|S><big|sum><rsub|s=1><rsup|S>h<around*|(|x<rsup|<around*|(|s|)>>|)>
    w<around*|(|x<rsup|<around*|(|s|)>>|)>>>>>
  </eqnarray>

  where <math|<around*|{|x<rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>> have
  been sampled from <math|q>. The function
  <math|w<around*|(|x|)>=<frac|p<around*|(|x|)>|q<around*|(|x|)>>> is
  referred to as the <with|font-shape|italic|weight function> in the
  literature. This estimator is unbiased, so that given enough samples we can
  accurately compute the expectation. The key feature that determines the
  efficiency of the method is the variance of the estimator. In the ideal
  case with <math|p=q> the weights will be one for every sample. This is the
  case were we can actualy sample from <math|p>, so we end up with <math|S>
  equally weighted samples. If the weights are highly variable however,
  certain samples will contribute significantly more than others to the
  estimator. As a general rule we want <math|q> to be a good match to
  <math|p> so the weights are nearly one.

  We have a great deal of flexbility in choosing <math|q>. The main
  constrains are

  <\itemize-dot>
    <item><math|q<around*|(|x|)>\<gtr\>0> whenever
    <math|p<around*|(|x|)>\<gtr\>0> so the weights are well defined

    <item>The support of <math|q> contains the support of <math|p>. In other
    words we can propose all possible values of <math|p> from <math|q>.
  </itemize-dot>

  <subsubsection|Sequential importance sampling>

  <\algorithm>
    <label|alg:sis><math|x<rsub|1>\<sim\>q<rsub|1><around*|(|\<cdot\>|)>>

    <math|w<rsub|1><around*|(|x<rsub|1>|)>\<leftarrow\><frac|\<mu\><around*|(|x<rsub|1>|)>
    g<around*|(|x<rsub|1>,y<rsub|1>|)>|q<rsub|1><around*|(|x<rsub|1>|)>>>

    for <math|t\<in\><around*|{|2,\<ldots\>,T|}>> do

    <\indent>
      <math|x<rsub|t>\<sim\>q<rsub|t><around*|(|\<cdot\>|)>>

      <math|w<rsub|t><around*|(|x<rsub|t>|)>\<leftarrow\><frac|f<around*|(|x<rsub|t-1>,x<rsub|t>|)>
      g<around*|(|x<rsub|t>,y<rsub|t>|)>|q<rsub|t><around*|(|x<rsub|t>|)>>>
    </indent>
  </algorithm>

  Identifying a good distribution to use for <math|q> when using IS can be
  hard if the dimensionality of our distribution is high. For example, if we
  imagine trying to use importance sampling for the hidden states of an HMM
  with length <math|T> we need to sample vectors of length <math|T>. Let
  <math|\<mu\><around*|(|x<rsub|1>|)>> be the probability of observing
  initial state <math|x<rsub|1>>, <math|f<around*|(|x<rsub|t>,x<rsub|t+1>|)>>
  be the probability of transitioning from state <math|x<rsub|t>> to
  <math|x<rsub|t+1>> and <math|g<around*|(|x<rsub|t>,y<rsub|t>|)>> be the
  probability of observing the value <math|y> given <math|x>. Then we have
  the following decomposition of the state probabilities.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<b-x\>|)>>|<cell|=>|<cell|\<mu\><around*|(|x<rsub|1>|)>g<around*|(|x<rsub|1>,y<rsub|1>|)>
    <big|prod><rsub|t=2><rsup|T>f<around*|(|x<rsub|t-1>,x<rsub|t>|)>
    g<around*|(|x<rsub|t>,y<rsub|t>|)>>>>>
  </eqnarray>

  Here we use the notation <math|\<b-x\>=<around*|(|x<rsub|1>,\<ldots\>,x<rsub|T>|)>>
  and in future we will also use <math|\<b-x\><rsub|u:v>=<around*|(|x<rsub|u>,x<rsub|u+1>,\<ldots\>,x<rsub|v-1>,x<rsub|v>|)>>
  so that <math|\<b-x\><rsub|1:t>=<around*|(|x<rsub|1>,\<ldots\>,x<rsub|t>|)>>
  and <math|\<b-x\>=\<b-x\><rsub|T>>. This suggests a sequential way to
  construct a proposal distribution. Rather try to propose <math|\<b-x\>>
  directly we will instead define a sequence of distributions
  <math|q<rsub|t>> to propose <math|x<rsub|t>>, so that we have\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|q<around*|(|\<b-x\>|)>>|<cell|=>|<cell|<big|prod><rsub|t=1><rsup|T>q<rsub|t><around*|(|x<rsub|t>|)>>>>>
  </eqnarray>

  \;

  With these definitions we can use Algorithm <reference|alg:sis> perform
  sequential importance sampling (SIS). This algorithm will produce the
  variable we wish to sample <math|\<b-x\>> and a collection of weights
  <math|\<b-w\>=<around*|(|w<rsub|1>,\<ldots\>,w<rsub|T>|)>>. The usefullness
  of this procedure stems from the following identity.

  <\eqnarray>
    <tformat|<table|<row|<cell|<big|prod><rsub|t=1><rsup|T>w<rsub|t><around*|(|x<rsub|t>|)>>|<cell|=>|<cell|<frac|\<mu\><around*|(|x<rsub|1>|)>
    g<around*|(|x<rsub|1>,y<rsub|1>|)>|q<rsub|1><around*|(|x<rsub|1>|)>><big|prod><rsub|t=2><rsup|T><frac|f<around*|(|x<rsub|t-1>,x<rsub|t>|)>
    g<around*|(|x<rsub|t>,y<rsub|t>|)>|q<rsub|t><around*|(|x<rsub|t>|)>>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<mu\><around*|(|x<rsub|1>|)>
    g<around*|(|x<rsub|1>,y<rsub|1>|)> <big|prod><rsub|t=2><rsup|T>f<around*|(|x<rsub|t-1>,x<rsub|t>|)>
    g<around*|(|x<rsub|t>,y<rsub|t>|)>|<big|prod><rsub|t=1><rsup|T>q<rsub|t><around*|(|x<rsub|t>|)>>>>|<row|<cell|>|<cell|=>|<cell|<frac|p<around*|(|\<b-x\>|)>|q<around*|(|\<b-x\>|)>>>>|<row|<cell|>|<cell|\<assign\>>|<cell|w<around*|(|\<b-x\>|)>>>>>
  </eqnarray>

  In other words we can compute the importance weight of <math|\<b-x\>>,
  <math|w<around*|(|\<b-x\>|)>>, by taking the product of the sequential
  importance weights <math|w<rsub|t><around*|(|x<rsub|t>|)>>. This allows us
  to use a simpler sequence of proposals <math|q<rsub|t>> on a lower
  dimensional space.\ 

  The proposal functions can depend on the previous values we have samples,
  so we have <math|q<rsub|t><around*|(|x<rsub|t>\|\<b-x\><rsub|1:t-1>|)>>. In
  the context of HMMs this means we could have
  <math|q<rsub|1><around*|(|x<rsub|1>|)>=\<mu\><around*|(|x<rsub|1>|)>> and
  <math|q<rsub|t><around*|(|x<rsub|t>\|\<b-x\><rsub|1:t-1>|)>=f<around*|(|x<rsub|t-1>,x<rsub|t>|)>>
  i.e. we propose from the initial state function and transition function.
  This commonly called the <with|font-shape|italic|bootstrap> proposal in the
  literature. A more efficient proposal is to use
  <math|q<rsub|t><around*|(|x<rsub|t>|)>=p<around*|(|x<rsub|t>\|\<b-x\><rsub|1:t-1>|)>=<frac|p<around*|(|\<b-x\><rsub|1:t>|)>|p<around*|(|\<b-x\><rsub|1:t-1>|)>>=<frac|p<around*|(|\<b-x\><rsub|1:t>|)>|<big|sum><rsub|x<rsub|t>>p<around*|(|\<b-x\><rsub|1:t>|)>>>
  which we refer to as the fully adapted proposal. The benefit of this
  proposal is tht it considers all the previous states and the data seen so
  far. The bootstrap proposal in contrast only considers the previous state,
  and ignores the data.

  <subsubsection|Sequential Monte Carlo>

  <\algorithm>
    <label|alg:smc>for <math|s\<in\><around*|{|1,\<ldots\>,S|}>> do

    <\indent>
      <math|<wide|x|~><rsub|1><rsup|<around*|(|s|)>>\<sim\>q<rsub|1><around*|(|\<cdot\>|)>>

      <math|w<rsub|1><rsup|<around*|(|s|)>>\<leftarrow\><frac|\<gamma\><rsub|1><around*|(|<wide|x|~><rsup|<around*|(|s|)>>|)>|q<rsub|1><around*|(|<wide|x|~><rsup|<around*|(|s|)>>|)>>>
    </indent>

    <math|<wide|w|\<bar\>><rsub|1><rsup|<around*|(|s|)>>\<leftarrow\><frac|w<rsub|1><rsup|<around*|(|s|)>>|<big|sum><rsub|s=1><rsup|S>w<rsub|1><rsup|<around*|(|s|)>>>>

    <math|<around*|{|x<rsub|1><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>\<leftarrow\><text|RESAMPLE><around*|(|<around*|{|<wide|x|~><rsub|1><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>,<around*|{|<wide|w|\<bar\>><rsub|1><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>|)>>

    for <math|t\<in\><around*|{|2,\<ldots\>,T|}>> do

    <\indent>
      for <math|s\<in\><around*|{|1,\<ldots\>,S|}>> do

      <\indent>
        <math|<wide|\<b-x\>|~><rsub|1:t><rsup|<around*|(|s|)>>\<sim\>q<rsub|t><around*|(|\<cdot\>\|\<b-x\><rsub|1:t-1><rsup|<around*|(|s|)>>|)>>

        <math|w<rsub|t><rsup|<around*|(|s|)>>\<leftarrow\><frac|\<gamma\><rsub|t><around*|(|<wide|\<b-x\>|~><rsub|1:t><rsup|<around*|(|s|)>>|)>|q<rsub|t><around*|(|<wide|\<b-x\>|~><rsub|1:t><rsup|<around*|(|s|)>>\|\<b-x\><rsub|1:t-1><rsup|<around*|(|s|)>>|)>>>
      </indent>

      <math|<wide|w|\<bar\>><rsub|t><rsup|<around*|(|s|)>>\<leftarrow\><frac|w<rsub|t><rsup|<around*|(|s|)>>|<big|sum><rsub|s=1><rsup|S>w<rsub|t><rsup|<around*|(|s|)>>>>

      <math|<around*|{|\<b-x\><rsub|1:t><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>\<leftarrow\><text|RESAMPLE><around*|(|<around*|{|<wide|\<b-x\>|~><rsub|1:t><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>,<around*|{|<wide|w|\<bar\>><rsub|t><rsup|<around*|(|s|)>>|}><rsub|s=1><rsup|S>|)>>
    </indent>
  </algorithm>

  The main problem with SIS is that the weights can be highly variable which
  leads to a poor estimate of the expectation. The basic issue is that each
  sample indepdently takes a path through the stat space and some will
  produce samples which are poor fit to <math|p>. The sequential Monte Carlo
  (SMC) algorithm outline in Algorithm <reference|alg:smc> seeks to address
  this problem by using a set of <math|S> interacting
  <with|font-shape|italic|particles>. Intuitively each particle is like a
  single sample from the SIS algorithm. The main difference is that we
  resample the particles at each step <math|t> with a probability
  proportional to their weights. The resampling step acts to prune aways
  particles with low weights and replace them by particles with higher
  weights. We can then estimate our expectation as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|p><around*|[|h|]>>|<cell|\<approx\>>|<cell|<big|sum><rsub|s=1><rsup|S>h<around*|(|\<b-x\><rsub|1:T><rsup|<around*|(|s|)>>|)>>>>>
  </eqnarray>

  where the weights do not appear because of the resampling step. If we
  neglect the final resampling step, then we would instead use the following
  estimator.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|p><around*|[|h|]>>|<cell|\<approx\>>|<cell|<big|sum><rsub|s=1><rsup|S>w<rsub|t><around*|(|<wide|\<b-x\>|~><rsub|1:T><rsup|<around*|(|s|)>>|)>h<around*|(|<wide|\<b-x\>|~><rsub|1:T><rsup|<around*|(|s|)>>|)>>>>>
  </eqnarray>

  There are different ways to resample the particles. The main constraint is
  that we must produce unbiased samples. Informally this means the
  probability of observing a particle after resampling should be proportional
  to its weight. The simplest scheme is multinomial sampling, were one draws
  the particle indices to keep from a multinomial distribution with
  parameters <math|S> and <math|<wide|\<b-w\>|\<bar\>><rsub|t>=<around*|(|<wide|w|\<bar\>><rsub|t><rsup|<around*|(|1|)>>,\<ldots\>,w<rsup|<around*|(|S|)>><rsub|t>|)>>.
  Other approaches such as stratified and systematic sampling can lead to
  lower variance estimates.

  Algorithm <reference|alg:smc> for SMC is given in a more general form than
  Algorithm <reference|alg:sis> for SIS. We no longer assume we have an HMM
  model. The assumption we know make is that the distribution were are
  interesed in takes the following form <math|p<around*|(|\<b-x\>|)>=<frac|\<pi\><around*|(|\<b-x\>|)>|Z>>.
  Here <math|\<pi\><around*|(|\<b-x\>|)>> is the unnormalized probability and
  <math|Z> the normalisation constant. In the Bayesian setup
  <math|\<pi\><around*|(|\<b-x\>|)>> would be the joint distribution. We then
  define a sequence of target distribution
  <math|<around*|{|\<gamma\><rsub|t><around*|(|\<b-x\><rsub|1:t>|)>|}><rsub|t=1><rsup|T>>
  as the numerator in our weight computation. The key constraint we have is
  that <math|\<gamma\><rsub|T><around*|(|\<b-x\><rsub|1:T>|)>\<propto\>p<around*|(|\<b-x\><rsub|1:T>|)>>,
  that is at the final step our target distribution has to be proportional to
  the distribution we are interested in sampling from.

  <subsubsection|Particle Metropolis-Hastings>

  <subsubsection|Particle Gibbs>

  <subsection|Hamiltonian Monte Carlo>

  The MH algorithm works by performing a biased random walk around the
  parameter space. The computational efficiency of this method depends
  critically on how well the chain explores the space. This is largely
  dictated by the proposal distribtion <math|q>. A seemingly strange feature
  of <math|q> is that it typically uses no information about the target
  distribution <math|p>. If we think in terms of optimization for a momemnt,
  then a random walk is not a very efficient way to find and optimum.
  Gradient based methods will typically be much better as they can leverage
  information about the direction of steepest ascent.

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
    <associate|alg:sis|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|alg:smc|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-1|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-2|<tuple|1.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-3|<tuple|1.1.1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-4|<tuple|1.1.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-5|<tuple|1.1.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-6|<tuple|1.1.4|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-7|<tuple|1.1.5|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-8|<tuple|1.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      1.<space|2spc>Miscallaneous material
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Sequential Monte Carlo
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|2tab>|1.1.1.<space|2spc>Important sampling
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|1.1.2.<space|2spc>Sequential importance
      sampling <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.1.3.<space|2spc>Sequential Monte Carlo
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>
    </associate>
  </collection>
</auxiliary>