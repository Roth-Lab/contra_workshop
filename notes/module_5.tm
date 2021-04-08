<TeXmacs|1.99.13>

<style|<tuple|tmbook|british|old-dots>>

<\body>
  <section|Advanced Inference Techniques>

  In this module we briefly review some more advanced techniques for fitting
  models. The main purpose is to highlight these approaches and their
  benefits, rather than provide a detailed description.

  <subsection|Variational inference>

  MCMC and MAP estimation represent two extremes of the possibilities for
  approximating a posterior distribution. MCMC is accurate but
  computationally expensive, whereas MAP is very inaccurate (in terms of
  approximating the whole distribution) but relatively fast. In between these
  extremes there are many other strategies for approximating posteriors.
  Variational inference (VI) has a emerged as one promising approach and has
  developed a fairly active research community <cite|blei2017variational>.

  VI frames the problem of approximating the posterior distribution as an
  optimisation problem. The basic idea is to identify a family of
  approximating distributions <math|q<around*|(|\<theta\>\|\<lambda\>|)>>
  which are in some sense simpler than the posterior
  <math|p<around*|(|\<theta\>\|X|)>>. The goal is then to find the settings
  of the variational parameters <math|\<lambda\>> that minimise the distance
  between <math|q> and the posterior. Typically the Kullback-Leibler
  divergence is used, though other alternatives are being explored.

  VI is typically much faster than MCMC methods and gives a more accurate
  posterior approximation than MAP estimation. Historically mean field VI
  (MFVI) has been the most popular approach <cite|beal2003variational>. MFVI
  assumes <math|q<around*|(|\<theta\>\|\<lambda\>|)>> decomposes as product
  of distribution for each parameter in the model i.e.
  <math|q<around*|(|\<theta\>\|\<lambda\>|)>=<big|prod><rsub|i=1><rsup|D>q<around*|(|\<theta\><rsub|i>\|\<lambda\><rsub|i>|)>>
  where <math|\<theta\>=<around*|(|\<theta\><rsub|1>,\<ldots\>,\<theta\><rsub|D>|)>>.
  The main issue with this approach is that it cannot capture correlation
  between variables. It was also hard to apply when the model was not
  composed of distributions in the conjugate exponential family. The primary
  reason for this constraint was the need to compute expectation under the
  approximating distribution.\ 

  Recent work has significantly relaxed the constraints of MFVI. The key
  insight has been that difficult to compute expectations can be approximated
  using Monte Carlo (MC) methods during the optimisation process.
  Specifically we use MC to estimate the gradient of the distance between
  <math|q<around*|(|\<theta\>\|\<lambda\>|)>> and
  <math|p<around*|(|\<theta\>\|X|)>> to perform
  <with|font-shape|italic|stochastic gradient descent>. This greatly expands
  the scope of models that VI can be applied to, and allows for approximating
  distributions which can capture correlation between model parameters. The
  key challenge is controlling the variance of the MC estimates. These
  advances have motivated an interesting area of research is using neural
  networks (deep learning) as part of the approximating distribution. Neural
  networks are powerful tools for function approximation and can in principle
  lead to quite accurate posterior approximation.

  Stochastic gradient descent has also been used to sub-sample data to
  estimate the gradient of the distance. This can lead to a significant
  reduction in computational complexity for each step of the gradient descent
  procedure. One example of the utility of this approach is fitting topic
  models to large data sets like Wikipedia with millions of data points.
  Sub-sampling turns out to be a major advantage for VI and is something that
  cannot be currently done with MCMC.\ 

  <subsection|Advanced MCMC methods>

  The Metropolis-Hastings (MH) algorithm is a widely used and often effective
  MCMC technique. However, it can struggle with high dimensional problems.
  The basic issue is the <with|font-shape|italic|curse of dimensionality>. In
  general if we propose new values far from the current one with MH they will
  be rejected. As a result we are forced to explore a small local space
  around our current value. When the dimensionality of the parameter space is
  high, this means that we cannot move any single component very far. One
  solution we have discussed is to <with|font-shape|italic|block> the data an
  perform updates on small subsets of variables. This approach begins to fail
  when there is strong correlation between variables in different blocks.

  The other common problem for MH is local optima in the posterior
  distribution. In general it can be hard for MH to move between modes in the
  distribution. This can be a problem if we initialise the MH sampler near a
  local mode which is far from the global optima. In this case we will get
  stuck sampling around the local mode and fail to explore the rest of the
  space where most distribution mass exists.

  There are a number of more advanced MCMC methods which can be used to
  overcome these issues. Many of these methods can be used together and with
  MH to produce better mixing samplers.

  <subsubsection|Simulated annealing>

  Simulated annealing (SA) is method from the optimisation literature
  inspired by a physical model of cooling metals. The basic idea of SA is to
  introduce a sequence of distributions. As we move through the sequence of
  distributions sampling from them becomes harder, until we reach the
  distribution we are interested in.\ 

  One way to construct this sequence is to anneal the likelihood as follows

  <\eqnarray>
    <tformat|<table|<row|<cell|p<rsub|\<beta\>><around*|(|\<theta\>\|X|)>>|<cell|=>|<cell|p<around*|(|\<theta\>\|X|)><rsup|\<beta\>>
    p<around*|(|\<theta\>|)><rsup|1-\<beta\>>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<frac|p<around*|(|X\|\<theta\>|)>|p<around*|(|X|)>>|)><rsup|\<beta\>>p<around*|(|\<theta\>|)>>>>>
  </eqnarray>

  where we take the sequence <math|\<beta\><rsub|1>=0\<less\>\<beta\><rsub|2>\<less\>\<cdots\>\<less\>\<beta\><rsub|T>=1>.
  The first distribution in the sequence is the prior and the last the
  posterior. To use this sequence in an MCMC framework we start at
  <math|\<beta\><rsub|1>> and run some form of MCMC update like MH. After
  some number of iterations we then move to <math|\<beta\><rsub|2>> keeping
  the last value from our previous updates. We continue this until we reach
  <math|\<beta\><rsub|T>> at which point we are sampling from the posterior.
  We then collect samples from this last iteration as our posterior
  approximation.

  The intuition behind SA is that it will flatten the posterior out early on,
  making it easy to move between modes. As the the annealing parameter
  increases, the distribution <math|p<rsub|\<beta\>>> will be become more
  rugged until we reach the posterior. Hopefully by that point the sampler is
  near the mode and sampling well from it. In principle SA can escape local
  optima from bad initialisation.\ 

  Tuning the values of <math|\<beta\><rsub|i>> and selecting the correct
  number is the main challenge for SA. If the distance between
  <math|\<beta\><rsub|i>> and <math|\<beta\><rsub|i+1>> is large, then the
  associated distributions will be quite different. This in turn can lead the
  algorithm to get trapped in modes.

  <subsubsection|Hamiltonian Monte Carlo>

  Hamiltonian Monte Carlo (HMC) is an efficient method for sampling high
  dimensional continuous parameters. Because HMC jointly updates many
  parameters (potentially all of them) it can work well when there is
  correlation. When the parameters are highly correlated HMC can
  significantly outperform MH.

  The basic idea is to propose a new value for the MCMC sampler by
  deterministically moving along a path defined by Hamiltonian dynamics. The
  connection to posterior inference is that the joint distribution appears as
  the potential energy term in the Hamiltonian. If we could follow the
  Hamiltonian path exactly we would then generate a new value which is a
  sample from our posterior distribution. In practice we cannot follow this
  path exactly, so we need to use numerical methods to approximate a set of
  differential equations (DEs). To correct for the numerical errors this
  introduces, we then need to perform an accept/reject step like MH.

  There are two parameters that need to be tuned: the step size in our
  approximation of the DE and the number of steps we use. Recent work
  developing automated approaches for these parameters has made HMC much
  easier to implement. One of the major applications has been in the
  probabilistic programming language STAN. Auto tuned HMC is the main
  algorithm used in STAN, and seems to work well for a wide range of models.

  The main deficiency of HMC is that it only applies to continuous
  parameters. There are a few ways around this issue. The first is to
  marginalise the discrete variables in the model. For example we can sum out
  the cluster indicator variables in mixture models, or use the
  forward-backward algorithm to marginalise the hidden states in an HMM. If
  this is not possible, HMC can still be used to update blocks of continuous
  variables in the same way as MH.

  <subsubsection|Sequential Monte Carlo>

  Sequential Monte Carlo (SMC) is a very general algorithm for sampling from
  high dimensional distributions. The basic idea is to break the sampling
  problem down into smaller and easier to sample problems. To achieve this
  SMC maintains a collection of <with|font-shape|italic|particles>
  representing a partial state. At each iteration of the SMC algorithm we
  update the particles and compute a weight for them. The weight, roughly
  speaking, is a measure of how good a fit the partial solution represented
  by the particle is to the data. Resampling is used after each iteration or
  when some criteria is met to replace low weight particles with high weight
  ones. At the end of the SMC algorithm we have a set of particles
  representing the sampled parameters, and a collection of weights that can
  be used to form an approximation to the posterior.

  SMC is often used to sample from distributions with an obvious sequential
  structure. The classic example is the posterior distribution of the hidden
  states of an HMM. In this case the <math|t<rsup|th>> iteration of the SMC
  algorithm proposes a new value for the hidden state <math|z<rsub|t>> of the
  chain. The weight function depends on the partial likelihood of the data up
  to <math|t>. While it is most obvious how to apply SMC in models with
  sequential structure, it is much more general. Interesting examples include
  inference in Bayesian mixture models, phylogenetics and probabilistic
  graphical models.

  There are a few considerations when designing an efficient SMC sampler. The
  first is to identify a sequential decomposition of the problem. In some
  cases like HMMs this is obvious, while others require more thought. The
  next consideration is the sequence of target distributions
  <math|<around*|{|\<gamma\><rsub|t>|}><rsub|t=1><rsup|T>>. This is related
  to sequential decomposition, but there is some flexibility. The only real
  constraint is that the final distribution in the sequence is proportional
  to the distribution we wish to sample i.e.
  <math|\<gamma\><rsub|T><around*|(|\<b-theta\><rsub|1:T>|)>\<propto\>p<around*|(|\<b-theta\><rsub|1:T>\|X|)>>
  for posterior inference. Here we use the notation
  <math|\<b-theta\><rsub|1:t>=<around*|(|\<theta\><rsub|1>,\<ldots\>,\<theta\><rsub|t>|)>>
  to indicate the sequence of parameters that have been sampled. The final
  design choice the set of proposal distributions
  <math|<around*|{|q<rsub|t>|}><rsub|t=1><rsup|T>> to use at each step. The
  only real constraint is that we must be able to propose every possible
  value of <math|\<theta\><rsub|t>> at each step.

  <subsubsection|Particle Gibbs>

  SMC can be a powerful method for sampling high dimensional parameters. It
  is especially useful when the parameters are discrete and HMC cannot be
  applied. One major deficiency of SMC is that we can only make a single
  pass. For example, if we want to update the hidden states of an HMM we have
  to fix the values of the other parameters like the transition matrix.\ 

  A simple idea is to use the blocking strategy like we do for MH. We could
  then use SMC to update the hidden states, and some other MCMC move to
  update the remaining parameters. We then iterate these steps to draw
  multiple samples from the posterior. Naively applying this approach does
  not lead to a valid sampler. The key problem is that SMC does not sample
  from the correct conditional distribution. The Particle Gibbs (PG)
  algorithm provides a solution to this problem.

  The key modification the PG algorithm makes to the proposed blocked
  sampler, is to use conditional SMC (cSMC) in place of SMC. This simply
  amounts to ensuring there is one particle which follows a path that leads
  to the current state. The modifications required to the standard SMC
  algorithm to achieve this are fairly minimal. Intuitively the conditional
  particle acts to keep the other particles in the SMC swarm near the current
  value. This is achieved through the resampling step, whereby the
  conditional particle can replace other particles.

  <\remark>
    Recall that in blocked MH sampling we have the parameter vector
    <math|\<theta\>=<around*|(|\<theta\><rsub|1>,\<theta\><rsub|2>|)>> and we
    would like to sample from the conditional distribution
    <math|p<around*|(|\<theta\><rsub|1>\|X,\<theta\><rsub|2>|)>> like a Gibbs
    sampler. In principle it seems like we need to know
    <math|p<around*|(|X,\<theta\><rsub|2>|)>> but this cancels away in the
    acceptance ratio so we only ever need to evaluate the joint distribution
    <math|p<around*|(|X,\<theta\><rsub|1>,\<theta\><rsub|2>|)>>. The PG
    algorithm essentially mimics this cancellation process through the
    addition of the conditional path.
  </remark>

  <\remark>
    The most challenging aspect of implementing the PG algorithm is typically
    doing the book keeping for the conditional path.
  </remark>

  <subsubsection|Parallel tempering>

  Parallel tempering (PT) is a very general approach to address the issues of
  multi-modal or hard to sample from posteriors. Like SA we have sequence of
  distributions <math|<around*|{|p<rsub|\<beta\><rsub|i>>|}><rsub|i=1><rsup|T>>.
  The difference is we concurrently perform updates for all these
  distributions using an MCMC algorithm. Periodically we propose to swap the
  parameters between two chains. The posterior approximation is generated by
  sampling from the chain with annealing parameter <math|\<beta\><rsub|T>=1>.

  The basic idea of PT is that chains with lower values of <math|\<beta\>>
  are able to explore the space more easily and jump between modes. By
  swapping values between chains we can slowly percolate new values to the
  top chain with parameter <math|\<beta\><rsub|T>> \ which approximates the
  posterior. This can lead to big jumps in the space which can cross between
  modes.

  The key issue is designing a good sequence of
  <math|\<beta\><rsub|1>\<less\>\<cdots\>\<less\>\<beta\><rsub|T>> and
  ensuring we have the right number of chains. If we do not use enough chains
  the probability of accepting a swap becomes low because the distributions
  are too different between chains. If we use too many chains, then it will
  take a very long time for new values to move from the lowest chain to the
  top. Historically this was a major challenge when implementing PT.\ 

  Recent work on <with|font-shape|italic|non-reversible> PT methods has made
  this less of an issue. The key insight from this work is that using a
  non-reversible scheme for swapping can drastically accelerate the
  probability of a parameter moving from the bottom chain to the top. In
  practice this means that we can use as many chains as we can afford
  computationally. An attractive feature of PT is that each chain can run
  independently on a separate CPU which allows for massive parallel
  computation. Only when swaps are proposed do we need to work serially.

  <subsubsection|Approximate Bayesian computation>

  Approximate Bayesian computation (ABC) is fairly recent approach for
  sampling from complex models where the likelihood may be intractable or
  expensive. These type of models come up frequently in population genetics.
  A simple example is the Wright-Fisher model. While it is reasonably easy to
  simulate from this model, it is hard to evaluate the probability of the
  observed data. This poses a problem for classic Bayesian inference where
  the likelihood is of central importance.

  The solution ABC applies for estimating model parameters in these models is
  surprisingly simple. Assume we have a model with parameters
  <math|\<theta\>>. We are able to simulate data given <math|\<theta\>>, but
  cannot easily evaluate the likelihood <math|p<around*|(|X\|\<theta\>|)>>.
  In ABC we simply propose values of <math|<wide|\<theta\>|~>> from the prior
  distribution, then simulate data <math|<wide|X|~>> and then compute a
  measure of distance between the simulated data <math|<wide|X|~>> and
  observed data <math|X>. If the distance is less than a threshold,
  <math|\<varepsilon\>>, we add the <math|<wide|\<theta\>|~>> to our
  collection of samples used to approximate the posterior.

  The accuracy of the posterior approximation depends on the value
  <math|\<varepsilon\>>. If <math|\<varepsilon\>> is large we will accept
  many samples, but have a less accurate posterior approximation. In contrast
  if we use a small value of <math|\<varepsilon\>> we will often reject
  samples and need to run more trials, but the approximation will be better.\ 

  One key design choice is the distance metric used to compare the observed
  data to the simulated data. If we have high dimensional data, then a metric
  which uses all the data will typically have large distances. Thus, in
  practice the data is often summarised using lower dimensional summary
  statistics which are then compared. For example, we could use the mean and
  variance of the data as our summary statistics and compute the euclidean
  distance. ABC is a fairly recent development, and there is a great deal of
  research going on in the field. One of the core problems being considered
  is how best to choose summary statistics.

  ABC is computationally demanding. Typically many millions of simulations
  will be required before a sufficient number of samples can be collected.
  Thus, while it could in principle be used to fit any model where other MCMC
  approaches apply, it will be far too computationally demanding.

  <\bibliography|bib|tm-plain|references>
    <\bib-list|2>
      <bibitem*|1><label|bib-beal2003variational>Matthew<nbsp>James Beal
      et<nbsp>al. <newblock><with|font-shape|italic|Variational algorithms
      for approximate Bayesian inference>. <newblock>University of London
      London, 2003.<newblock>

      <bibitem*|2><label|bib-blei2017variational>David<nbsp>M Blei, Alp
      Kucukelbir<localize|, and >Jon<nbsp>D McAuliffe. <newblock>Variational
      inference: a review for statisticians.
      <newblock><with|font-shape|italic|Journal of the American Statistical
      Association>, 112(518):859\U877, 2017.<newblock>
    </bib-list>
  </bibliography>
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
    <associate|auto-10|<tuple|1.2.6|5>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|1>>
    <associate|auto-4|<tuple|1.2.1|2>>
    <associate|auto-5|<tuple|1.2.2|2>>
    <associate|auto-6|<tuple|1.2.3|2>>
    <associate|auto-7|<tuple|1.2.4|3>>
    <associate|auto-8|<tuple|1.2.5|3>>
    <associate|auto-9|<tuple|1.2.6|4>>
    <associate|bib-beal2003variational|<tuple|1|5>>
    <associate|bib-blei2017variational|<tuple|2|5>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      blei2017variational

      beal2003variational
    </associate>
    <\associate|toc>
      1.<space|2spc>Advanced Inference Techniques
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Variational inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2.<space|2spc>Advanced MCMC methods
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|1.2.1.<space|2spc>Simulated annealing
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.2.2.<space|2spc>Hamiltonian Monte Carlo
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.2.3.<space|2spc>Sequential Monte Carlo
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|2tab>|1.2.4.<space|2spc>Particle Gibbs
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|2tab>|1.2.5.<space|2spc>Parallel tempering
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|1.2.6.<space|2spc>Approximate Bayesian
      computation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-shape|<quote|small-caps>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <pageref|auto-10><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>