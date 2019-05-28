<TeXmacs|1.99.8>

<style|tmbook>

<\body>
  Mixture models are popular and useful class of probabilistic models. Put
  simply mixture models assume there is some unknown group structure among
  the data. The group a data point belongs to influences the observed value.

  To understand the basic problem imagine we have a two boxes of coin each
  with 100 coins. One box of coins is fair, that is the probability of
  observing heads or tails is 0.5 for either case. The other box of coins is
  biased so that the probability of observing heads is 0.4 and the
  probability of observing tails is 0.6. Now the coins are mixed up together
  in a single box.\ 

  <\problem*>
    You are given a random coin from the box and allowed to flip it 30 times.
    Can you predict whether the coin is biased?
  </problem*>

  <\solution>
    We will develop a simple probabilistic model to answer this questions.
    First we note that a natural way to model observing <math|x> heads from
    <math|n> coins flips when the probability of heads is <math|p> is to use
    a <math|<text|Binomial><around*|(|x\|n,p|)>> distribution. We also
    observe that there are 100 coins from each box, so the probability a
    random coin is from the fair coin box is
    <math|<frac|100|200>=<frac|1|2>=0.5> and from the biased box is also
    <math|0.5>.\ 

    Now we introduce some notation. Let

    <\itemize-dot>
      <item><math|X> be the number of coins which show heads

      <item><math|n> be the number of coin flips

      <item><math|p<rsub|0>=0.5> be the probability a fair coin shows heads

      <item><math|p<rsub|1>=0.4> be the probability a biased coin shows heads

      <item><math|Z> be a variable which equals 1 when the coins is biased
      and 0 otherwise
    </itemize-dot>

    Then we have the following hierachichal model

    <\eqnarray>
      <tformat|<table|<row|<cell|Z>|<cell|\<sim\>>|<cell|<text|Bernoulli><around*|(|\<cdot\>\|0.5|)>>>|<row|<cell|X\|n,p<rsub|0>,p<rsub|1>,Z=z>|<cell|\<sim\>>|<cell|<choice|<tformat|<table|<row|<cell|<text|Binomial><around*|(|\<cdot\>\|n,p<rsub|0>|)>>|<cell|<text|if>>|<cell|z=0>>|<row|<cell|<text|Binomial><around*|(|\<cdot\>\|n,p<rsub|1>|)>>|<cell|<text|if>>|<cell|z=1>>>>>>>|<row|<cell|>|<cell|=>|<cell|<text|><text|Binomial><around*|(|\<cdot\>\|n,p<rsub|z>|)>>>>>
    </eqnarray>

    Now the probability the coin is biased is

    <\eqnarray>
      <tformat|<table|<row|<cell|P<around*|(|Z=1\|X=x|)>>|<cell|=>|<cell|<frac|P<around*|(|X=x\|Z=1|)>
      P<around*|(|Z=1|)>|P<around*|(|X=x|)>>>>|<row|<cell|>|<cell|=>|<cell|<frac|P<around*|(|X=x\|Z=1|)>
      P<around*|(|Z=1|)>|P<around*|(|X=x\|Z=0|)>P<around*|(|Z=0|)>+P<around*|(|X=x\|Z=1|)>
      P<around*|(|Z=1|)>>>>>>
    </eqnarray>

    where the first line follows from Bayes' rule and the second from the law
    of total probability. From the definition of a Bernoulli we have that\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|P<around*|(|Z=1|)>>|<cell|=>|<cell|0.5>>|<row|<cell|P<around*|(|Z=0|)>>|<cell|=>|<cell|1-P<around*|(|Z=1|)>>>|<row|<cell|>|<cell|=>|<cell|0.5>>>>
    </eqnarray>

    And from out model assumptions we have

    <\eqnarray>
      <tformat|<table|<row|<cell|P<around*|(|X=x\|Z=1|)>>|<cell|=>|<cell|<binom|n|x>0.4<rsup|x>
      <around*|(|1-0.4|)><rsup|n-x>>>|<row|<cell|P<around*|(|X=x\|Z=1|)>>|<cell|=>|<cell|<binom|n|x>0.5<rsup|x>
      <around*|(|1-0.5|)><rsup|n-x>>>>>
    </eqnarray>

    Substituting these in we have\ 

    <\eqnarray>
      <tformat|<table|<row|<cell|P<around*|(|Z=1\|X=x|)>>|<cell|=>|<cell|<frac|<binom|n|x>0.4<rsup|x>
      <around*|(|1-0.4|)><rsup|n-x>\<times\>0.5|<binom|n|x>0.5<rsup|x>
      <around*|(|1-0.5|)><rsup|n-x>\<times\>0.5+<binom|n|x>0.4<rsup|x>
      <around*|(|1-0.4|)><rsup|n-x>\<times\>0.5>>>|<row|<cell|>|<cell|=>|<cell|<frac|0.4<rsup|x>
      <around*|(|1-0.4|)><rsup|n-x>|0.5<rsup|x>
      <around*|(|1-0.5|)><rsup|n-x>+0.4<rsup|x>
      <around*|(|1-0.4|)><rsup|n-x>\<times\>0.5>>>>>
    </eqnarray>
  </solution>

  <\exercise>
    What are the probabilities the coin is biased if

    <\enumerate-alpha>
      <item>We observe 10 heads?

      <item>We observe 15 heads?

      <item>We observe 25 heads?
    </enumerate-alpha>
  </exercise>

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