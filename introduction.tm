<TeXmacs|1.99.8>

<style|<tuple|generic|british>>

<\body>
  <section|Introduction>

  <verbatim|>The development of high throughput assays such as short read
  sequencing are driving cancer biology to become a data rich scientific
  field. Analysing this data and extracting biologically meaning is becoming
  increasingly challenging as large data sets filled with biological and
  technical noise are created. There are a many approaches computational
  biologists can use to meet this challenge. In these notes we will consider
  how the framework of probabilistic modelling can be useful. We will see
  that probabilistic modelling provides a coherent framework to analyse noisy
  data and extract interpretable results. We adopt the Bayesian view of
  statistics in these notes, where model parameters are treated as random
  variables. The Bayesian view provides a systematic framework to construct
  and fit probabilistic models.

  <subsection|Why do we need probabilistic models>

  Real world data is not perfect, there is always \Pnoise\Q whether it be
  from measurement error or some underlying variation. This is especially
  true in biology with the development of high throughput assays such as
  short read sequencing. As we develop new methods such as single cell
  sequencing there is a constant tension between throughput and accuracy. For
  example, when we sequence single cells we often have a choice. Do we
  sequence more cells with low coverage or fewer cells with high coverage.
  Lower coverage data will have more technical noise, but allow us to survey
  a broader population and account for biological variability. How can we
  make use of such data in the presence of this noise?

  In the simplest sense probabilistic modelling is about looking at the world
  when there is noise. More accurately it is about constructing models which
  account for some form of randomness. A common example across many
  scientific fields is to assume the observed values are affected by
  \Qrandom\Q noise. This assumption stems from issues such as measurement
  error. While we believe there is some true value we wish to measure, the
  instruments available to make measurements intrinsically perturb this
  value. We can't learn a great deal from a single measurement, but if we
  make a lot measurements then we can average them to get a more accurate
  estimate. Of course the confidence we have in this measurement will be
  dependent on the variation we observe in the data. Classically we report
  such uncertainty with quantities like confidence intervals. If we adopt
  such a process we are actually constructing a simple probabilistic model.
  We are assuming the data is sampled from a Normal distribution with a mean
  and variance parameter that we don't know. Roughly speaking our average is
  the estimate of the mean (the true parameter) and the confidence interval
  an estimate of the variance (noise).

  Probabilistic modelling formalises the intuition of the last paragraph. We
  construct models where we encode quantities of interest as parameters and
  observations as random samples from a distribution that depends on these
  parameters. Models can be complex and have 1000s of parameters and
  1,000,000s of observation. The challenge for computational biologists is
  how to construct and fit these models. In these notes we will look at few
  key problems.

  <\enumerate-numeric>
    <item>How do we construct an appropriate model for the data?

    <item>How can estimate the model parameters?

    <item>How do report these model parameters i.e. how certain are we?
  </enumerate-numeric>

  In many ways, the first point is the most challenging. Constructing useful
  probabilistic models requires an understanding of the problem domain. The
  most critical thing that needs to be identified is: <strong|What is the
  question?> While this may seem simple, precisely defining what you are
  trying to ask is the single most important step in data analysis and is
  often poorly done. Once you have a question identified you then need to
  leverage knowledge of the domain to begin constructing models. What are
  biologically realistic assumptions? What is the data and what are the types
  of noise it exhibits? Building models which encapsulate this information is
  as much an art as science. Models that are too simple risk fitting the data
  poorly, whereas models that are to complex risk becoming impossible to
  interpret. In the next section we will discuss some basic ingredients to
  building models.

  The second and third point are somewhat more straightforward to address. We
  will adopt a Bayesian viewpoint in these notes. The clear advantage of the
  Bayesian method is that it provides a simple prescription for how to
  estimate model parameters. <em|Everything you want to know is encapsulated
  in the posterior distribution>. The hard part then is not how to estimate
  the model parameters or report uncertainty, but how to compute the
  posterior. For reason to be discussed later this is computationally
  challenging.

  One factor that is often overlooked is the interplay between model
  construction and inference. If we build a model that we cannot fit, then it
  is useless. Conversely, if we do not know how to fit a model we will not
  construct it. Thus, having an extensive toolbox of methods for fitting
  models opens up broader modelling possibilities.

  <subsection|Probability distributions>

  The basic ingredient of probabilistic modelling is probability theory. As
  we will see later we will make assumptions about the distributions that
  govern the data and model parameters. To be able to do this we need to have
  a reasonable knowledge of the distributions out there. In this section we
  will review a few common distributions, though this is by no means an
  exhaustive. Before we do that we discuss a basic rule of thumb that can
  help guide the selection of distribution.

  <subsubsection|Picking a distribution>

  The most fundamental question that needs to be addressed when choosing a
  distribution for modelling is: <em|What possible values can the quantity
  have?> For example, if your observed data can only take values in the
  positive integers (such as read depth), then a Normal distribution will be
  a poor choice to use in the model. Recall the set of possible values
  (formally the support of the distribution) of a Normal random variable is
  <math|\<bbb-R\>>, all possible real numbers. A much better choice in this
  case is the Poisson distribution which can take values in
  <math|\<bbb-Z\><rsup|+>>, the set of all non-negative integers.

  This simple rule of thumb is surprisingly powerful, it at least guarantees
  the model makes some sense. Of course it is just a rule of thumb. There are
  many distributions which have support in <math|\<bbb-Z\><rsup|+>>, the
  Negative-Binomial (NB) being another popular choice. Why you would prefer
  the Poisson or NB then comes down to several other factors. The most
  important is how variable is the data? The NB is a two parameter
  distribution, whereas the Poisson has only a single parameter. The extra
  parameter in the NB provides more flexibility in modelling the variance of
  the data. In fact, the NB is what is called over-dispersed relative to the
  Poisson. There are many pairs of distributions which have this
  relationship. Perhaps the most famous pair is the Normal and Student-t
  (over-dispersed) distributions.\ 

  There are other factors such distribution skew i.e. is the distribution
  symmetric and tail behaviour i.e. do we have more outliers than we expect.
  Knowing which distribution to use is often an empirical question, and
  rigorously should be performed as part of model selection. In practice it
  is best to start with the simplest (fewest parameter) distribution, that is
  easiest to interpret. If you observe a poor fit to the data, for example
  your model predictions deviate wildly from expectation, then begin to
  consider more complex choices.

  <subsection|Bayesian inference>

  Once you have constructed a model you need to <em|fit> the model to the
  data. What fit means depends on your chosen paradigm. In statistics there
  are two major paradigms, the frequentist and Bayesian paradigms. There are
  philosophical and practical differences between then two. We adopt the
  Bayesian paradigm in this work. The main reasons for doing do are:

  <\itemize-dot>
    <item>Prior information - There is a natural mechanism to incorporate
    existing knowledge into the Bayesian framework through the use of prior
    distribution

    <item>Hierarchical models - Model parameters are random variables in the
    Bayesian paradigm. This means that they can have distributions which
    govern their behaviour. The parameters in these distributions can in turn
    have their own distributions and so on. This allows for the construction
    of powerful hierarchies of variables that can be exploited to <em|share
    statistical> strength.

    <item>Coherent inference - As we will see fitting model parameters in the
    Bayesian paradigm amounts to computing the posterior distribution.
  </itemize-dot>

  <subsubsection|The posterior distribution>

  Abstractly we assume that we have a set of observed data <math|X> and model
  with a collection of parameters <math|\<theta\>>. The observed data
  <math|X> is modelled by a distribution which depends on <math|\<theta\>>,
  often called the <em|likelihood>. We denote the density of this
  distribution <math|p<around*|(|X\|\<theta\>|)>>. In the Bayesian paradigm
  we also have a distribution for the parameters, referred to as the
  <em|prior distribution> <math|p<around*|(|\<theta\>|)>>. The prior
  distribution reflects our belief about the parameters before observing the
  data. Our goal is then to compute the <em|posterior distribution>
  <math|p<around*|(|\<theta\>\|X|)>>. The posterior distribution describes
  our belief about the parameter values after observing the data. Computing
  the posterior is a trivial exercise in applying Bayes' rule.

  <\eqnarray>
    <tformat|<table|<row|<cell|p<around*|(|\<theta\>\|X|)>>|<cell|=>|<cell|<frac|p<around*|(|X\|\<theta\>|)>
    p<around*|(|\<theta\>|)>|p<around*|(|X|)>><eq-number><label|eq:posterior>>>|<row|<cell|p<around*|(|X|)>>|<cell|=>|<cell|<big|int>p<around*|(|X\|\<theta\>|)>
    p<around*|(|\<theta\>|)> \<mathd\>\<theta\><eq-number><label|eq:norm>>>>>
  </eqnarray>

  While it is conceptually simple to apply equation <reference|eq:posterior>
  to obtain the posterior, it is computationally challenging. The main issue
  is the need to compute the integral in <reference|eq:norm> which is
  typically high-dimensional and has no closed form. We will see later how we
  can address this problem using maximum a posterior (MAP) or Monte Carlo
  Markov Chain (MCMC) methods to approximate the posterior.

  <\example>
    \;

    We consider the simple problem of inferring the probability a coin will
    show a head when flipped. Assume we have coin with probability
    <math|\<rho\>> of showing heads and we perform <math|n> coin flips,
    <math|x> of which show heads. A reasonable likelihood for this problem is
    a Binomial with parameters <math|n> and <math|\<rho\>>.\ 

    If we assume no existing knowledge then a reasonable prior distribution
    is a continuous distribution. More generally we can assume that we
    <math|\<rho\>> follows a Beta distribution with parameters <math|a> and
    <math|b>. If we set <math|a=b=1> we recover the Uniform distribution.
    Alternatively a value <math|a=b=2> would reflect a weak prior belief that
    <math|\<rho\>> is near 0.5 and the coin is fair.

    Thus the model we have is\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|\<rho\>>|<cell|\<sim\>>|<cell|<text|Beta><around*|(|\<cdot\>\|a,b|)>>>|<row|<cell|x\|\<rho\>>|<cell|\<sim\>>|<cell|<text|Binomial><around*|(|\<cdot\>\|n,\<rho\>|)>>>>>
    </eqnarray>

    if this notation is unfamiliar see the section on hierarchical modelling.
    Also note we are suppressing the dependency on <math|n> which observed
    and not random. We can write down the distributions\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<rho\>|)>>|<cell|=>|<cell|<frac|1|\<cal-B\><around*|(|a,b|)>>
      \<rho\><rsup|a-1> <around*|(|1-\<rho\>|)><rsup|b-1>
      \<bbb-I\><around*|(|\<rho\>\<in\><around*|[|0,1|]>|)>>>|<row|<cell|p<around*|(|x\|\<rho\>|)>>|<cell|=>|<cell|<binom|n|x>
      \<rho\><rsup|x> <around*|(|1-\<rho\>|)><rsup|n-x>>>>>
    </eqnarray>

    and we have the joint distribution

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|x,\<rho\>|)>>|<cell|=>|<cell|p<around*|(|x\|\<rho\>|)>
      p<around*|(|\<rho\>|)>>>|<row|<cell|>|<cell|=>|<cell|<binom|n|x>
      \<rho\><rsup|x> <around*|(|1-\<rho\>|)><rsup|n-x>
      <frac|1|\<cal-B\><around*|(|a,b|)>> \<rho\><rsup|a-1>
      <around*|(|1-\<rho\>|)><rsup|b-1>>>|<row|<cell|>|<cell|=>|<cell|<frac|<binom|n|x>
      |\<cal-B\><around*|(|a,b|)>>\<rho\><rsup|x+a-1>
      <around*|(|1-\<rho\>|)><rsup|n-x+b-1>>>>>
    </eqnarray>

    To compute the posterior we need to evaluate

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<rho\>\|x|)>>|<cell|=>|<cell|<frac|p<around*|(|x,\<rho\>|)>|<big|int>p<around*|(|x,\<rho\>|)>
      \<mathd\>\<rho\>>>>|<row|<cell|>|<cell|=>|<cell|<frac|<frac|<binom|n|x>
      |\<cal-B\><around*|(|a,b|)>>\<rho\><rsup|x+a-1>
      <around*|(|1-\<rho\>|)><rsup|n-x+b-1>|<big|int><frac|<binom|n|x>
      |\<cal-B\><around*|(|a,b|)>>\<rho\><rsup|x+a-1>
      <around*|(|1-\<rho\>|)><rsup|n-x+b-1>
      \<mathd\>\<rho\>>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<rho\><rsup|x+a-1>
      <around*|(|1-\<rho\>|)><rsup|n-x+b-1>|<big|int>\<rho\><rsup|x+a-1>
      <around*|(|1-\<rho\>|)><rsup|n-x+b-1> \<mathd\>\<rho\>>>>>>
    </eqnarray>

    Now the challenge is to evaluate the integral in the denominator. We note
    a trick here which that the term <math|\<rho\><rsup|x+a-1>
    <around*|(|1-\<rho\>|)><rsup|n-x+b-1>> is exactly the term that appears
    in a <math|<text|Beta><around*|(|x+a,n-x+b|)>> distribution neglecting
    the normalisation constant. So it follows the integral equals the
    normalisation constant <math|\<cal-B\><around*|(|x+a,n-x+b|)>>. So we
    have\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<rho\>\|x|)>>|<cell|=>|<cell|<frac|1|\<cal-B\><around*|(|x+a,n-x+b|)>>
      \<rho\><rsup|x+a-1> <around*|(|1-\<rho\>|)><rsup|n-x+b-1>>>>>
    </eqnarray>

    or alternatively

    <\eqnarray>
      <tformat|<table|<row|<cell|\<rho\>\|x>|<cell|\<sim\>>|<cell|<text|Beta><around*|(|\<cdot\>\|x+a,n-x+a|)>>>>>
    </eqnarray>
  </example>

  In the previous example we could explicitly compute the normalisation
  constant. This is possible because the Beta and Binomial distributions are
  a conjugate pair. Conjugacy plays an important in many Bayesian models,
  making exact computation of the posterior possible. However, there are many
  models which are non-conjugate and for which the normalisation cannot be
  compute explicitly. In the previous example this is still not a major
  problem since we could simply numerically compute the one dimensional
  integral over <math|\<rho\>>. This can be easily done with functions in
  your favourite scientific computing package. However, if the dimensionality
  of the integral is even moderately high this becomes computationally
  prohibitive and alternative approaches will be required.

  <subsubsection|Summarising posteriors>

  Once we have a posterior distribution (or approximation) the question
  becomes how to report the results. While the posterior distribution
  provides all the information we could want, it is difficult for humans to
  make sense of complex and multi-dimensional distributions. As a result we
  often resort to reporting some form of point or interval estimate to
  summarise the posterior. The simplest quantities to report are the
  posterior mean and variance. These are often reasonable values, but care
  must be taken if the posterior is multi-modal. In that case the mean may
  not be a typical value for the distribution. An alternative approach is to
  report an <math|\<alpha\>> credible regions. This is an region which
  contains <math|\<alpha\>> of the mass of the posterior. There are many ways
  to construct such intervals, one sensible way is to find the smallest
  region <math|R> where <math|<big|int><rsub|R>p<around*|(|\<theta\>\|X|)>\<mathd\>\<theta\>=\<alpha\>>.
  This is referred to as the high probability density (HPD) interval, since
  by definition it captures the regions of highest density.

  <\example>
    \;

    We return to the coin flipping example. From the previous section we now
    the posterior distribution is\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<rho\>\|x|)>>|<cell|=>|<cell|<frac|1|\<cal-B\><around*|(|x+a,n-x+b|)>>
      \<rho\><rsup|x+a-1> <around*|(|1-\<rho\>|)><rsup|n-x+b-1>>>>>
    </eqnarray>

    it is easy to check the mean of this distribution is
    <math|<frac|x+a|n+a+b>> and the variance is
    <math|<frac|x+a|<around*|(|n+a+b|)><rsup|2>>>.
  </example>

  A more general theory of Bayesian point estimation exists and makes use of
  the concept of a loss function. A loss function, <math|L<around*|(|x,y|)>>,
  is bi-variate function which takes on positive real value. The
  interpretation is that the functions represents our perceived loss if we
  were to estimate <math|x> but the true value is <math|y>. In a Bayesian
  setting we seek to minimise the expected loss under the posterior as our
  best estimate.\ 

  <\eqnarray>
    <tformat|<table|<row|<cell|<wide|\<theta\>|^>>|<cell|=>|<cell|<below|<text|argmin>|\<theta\><rprime|'>><above||>
    <big|int>L<around*|(|\<theta\><rprime|'>,\<theta\>|)>
    p<around*|(|\<theta\>\|X|)> \<mathd\>\<theta\>>>>>
  </eqnarray>

  Common examples of loss functions for continuous variables include the
  <math|L<rsup|1>> and <math|L<rsup|2>> norm

  <\eqnarray>
    <tformat|<table|<row|<cell|L<around*|(|x,y|)>>|<cell|=>|<cell|<around*|\||x-y|\|>>>|<row|<cell|L<around*|(|x,y|)>>|<cell|=>|<cell|<around*|\<\|\|\>|x-y|\<\|\|\>><rsup|2>>>>>
  </eqnarray>

  Other loss functions can be used to estimate summary regions or even
  distributions. The same basic framework also applies to discrete variables,
  though the minimisation of the loss can be more difficult if the space is
  large as we cannot use gradient methods.

  <subsubsection|Hierarchical models>

  The theory of Bayesian inference is coherent and elegantly contained in a
  few basic principles. However, the real utility of the paradigm is the
  ability to construct complex hierarchical models. This is possible because
  we view model parameters are random variables, in the same way we view
  data. There are several useful notational constructs for describing
  hierarchical models. The most explicit is to list out the distributional
  assumptions. For example, we have many statements of the form
  <math|a\|b\<sim\><text|G><around*|(|\<cdot\>\|b|)>>. This reads as follows:
  the distribution of variable <math|a> given variable <math|b> follows
  distribution <math|<text|G>> with parameter <math|b>. You may also see this
  written as <math|a\|b\<sim\><text|G><around*|(|a\|b|)>> or
  <math|a\|b\<sim\><text|G><around*|(|b|)>>. In general the distribution
  <math|<text|G>> can depend on <math|b> in more complex ways, possibly
  through some transformation as the next example illustrates.

  <\example>
    \;

    Here we will describe the basic Bayesian linear regression model. Recall
    that linear regression assumes we observe covariate data
    <math|\<b-X\>=<around*|(|\<b-x\><rsub|1>,\<ldots\>,\<b-x\><rsub|N>|)>>
    and outcome variables <math|\<b-y\>=<around*|(|y<rsub|1>,\<ldots\>,y<rsub|N>|)>>
    for <math|N> data points. Here <math|\<b-x\><rsub|n>\<in\>\<bbb-R\><rsup|D>>
    and <math|y<rsub|n>\<in\><around*|{|0,1|}>>. Informally we observe some
    real valued covariates for a data point <math|\<b-x\><rsub|n>>, and we
    believe these influence the outcome variables <math|y<rsub|n>> which
    takes on two values. We will also assume there is regression coefficient
    vector <math|\<b-beta\>=<around*|(|\<beta\><rsub|1>,\<ldots\>,\<beta\><rsub|d>|)>\<in\>\<bbb-R\><rsup|D>>
    that controls the influence of the carats on the outcome. A common use
    for this model is to learn predictors to separate data into two classes,
    such as patients that will respond to treatment (<math|y<rsub|n>=1>) and
    patients that will not (<math|y<rsub|n>=0>).\ 

    Since the variables <math|y<rsub|n>> are 0/1 valued a natural
    distribution to model them is the Bernoulli. We will assume the
    regression covariates are independent and model <math|\<beta\><rsub|d>>
    as Normally distributed. To link these two items together we need make
    the assumption the probability <math|y<rsub|n>=1> is of the form
    <math|\<sigma\><around*|(|<big|sum><rsub|d>\<beta\><rsub|d> x<rsub|n
    d>|)>> where <math|\<sigma\><around*|(|x|)>=<frac|1|1+e<rsup|-x>>> is the
    logistic function.\ 

    Our model then has the form

    <\eqnarray>
      <tformat|<table|<row|<cell|\<beta\><rsub|d>\|\<mu\>,\<sigma\><rsup|2>>|<cell|\<sim\>>|<cell|<text|Normal><around*|(|\<cdot\>\|\<mu\>,\<sigma\><rsup|2>|)>>>|<row|<cell|y<rsub|n>\|\<b-x\><rsub|n>,\<b-beta\>>|<cell|\<sim\>>|<cell|<text|Bernoulli><around*|(|\<cdot\>\|\<sigma\><around*|(|<big|sum><rsub|d>\<beta\><rsub|d>
      x<rsub|n d>|)>|)>>>>>
    </eqnarray>

    where <math|\<mu\>> and <math|\<sigma\><rsup|2>> are hype-parameters we
    assume known.
  </example>

  By listing out the distributional assumptions we fully specify the model.
  This notation also lays bare any conditional independence assumptions.
  Recall variables are conditionally independent if the value of one variable
  is independent from another given the value of a third. Conditionally
  independence is useful as it allows to construct modular models, where we
  need only focus on the distributional assumptions of a small number of
  variables. It can also have computational implications, for example
  allowing for the design of Gibbs samplers. Once we are given the
  conditional distributions, it is straightforward to write down the joint
  distribution for the model.

  <\example>
    \;

    Following the previous example we wish to write down the joint
    distribution

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<b-X\>,\<b-y\>,\<b-beta\>,\<mu\>,\<sigma\><rsup|2>|)>>|<cell|=>|<cell|<big|prod><rsub|n=1><rsup|N>p<around*|(|y<rsub|n>,\<b-x\><rsub|n>,\<b-beta\>,\<mu\>,\<sigma\><rsup|2>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|prod><rsub|n=1><rsup|N>p<around*|(|y<rsub|n>\|\<b-x\><rsub|n>,\<b-beta\>|)>
      p<around*|(|\<beta\><rsub|d>\|\<mu\>,\<sigma\><rsup|2>|)>>>|<row|<cell|>|<cell|=>|<cell|p<around*|(|\<beta\><rsub|d>\|\<mu\>,\<sigma\><rsup|2>|)><big|prod><rsub|n=1><rsup|N>p<around*|(|y<rsub|n>\|\<b-x\><rsub|n>,\<b-beta\>|)>>>|<row|<cell|>|<cell|=>|<cell|\<cal-N\><around*|(|\<beta\><rsub|d>\|\<mu\>,\<sigma\><rsup|2>|)><big|prod><rsub|n=1><rsup|N><around*|[|\<sigma\><around*|(|<big|sum><rsub|d>\<beta\><rsub|d>
      x<rsub|n d>|)>|]><rsup|y<rsub|n>><around*|[|1-\<sigma\><around*|(|<big|sum><rsub|d>\<beta\><rsub|d>
      x<rsub|n d>|)>|]><rsup|1-y<rsub|n>>>>>>
    </eqnarray>

    where <math|\<cal-N\>> indicates the Normal distribution density. Now to
    compute the posterior for <math|\<b-beta\>> we have

    <\eqnarray>
      <tformat|<table|<row|<cell|p<around*|(|\<b-beta\>\|\<b-X\>,\<b-y\>,\<mu\>,\<sigma\><rsup|2>|)>>|<cell|=>|<cell|<frac|p<around*|(|\<b-X\>,\<b-y\>,\<b-beta\>,\<mu\>,\<sigma\><rsup|2>|)>|<big|int>p<around*|(|\<b-X\>,\<b-y\>,\<b-beta\>,\<mu\>,\<sigma\><rsup|2>|)>
      \<mathd\>\<b-beta\>>>>>>
    </eqnarray>

    However, the integral in the denominator is no longer simple. If <math|D>
    is small then we could attempt to compute it numerically, but for
    realistic problems we will need to turn to alternative methods.
  </example>

  One other benefit of listing out the distributional assumptions is that it
  makes it possible to simulate data from the model. This can be particularly
  useful for exploring the data generated from model and quickly evaluating
  whether it looks to match the observed values. It can also provide a means
  to explore the behaviour of the model in different parameter regimes, such
  as high variance or varying numbers of data points. Finally, it can also be
  a powerful debugging tool to test inference algorithms such as MCMC
  methods. If the simulated data comes from a distribution that is not very
  noisy we should expect our inference algorithm to infer values that are
  close to the known truth.

  \;

  \ 

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
    <associate|auto-2|<tuple|1.1|1|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-3|<tuple|1.2|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-4|<tuple|1.2.1|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-5|<tuple|1.3|2|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-6|<tuple|1.3.1|3|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-7|<tuple|1.3.2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|auto-8|<tuple|1.3.3|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|eq:norm|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
    <associate|eq:posterior|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_11.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Why do we need probabilistic
      models <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Probability distributions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|2tab>|1.2.1<space|2spc>Picking a distribution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>Bayesian inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.3.1<space|2spc>The posterior distribution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>