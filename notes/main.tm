<TeXmacs|1.99.8>

<style|generic>

<\body>
  <include|module_2.tm>
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-type|letter>
    <associate|page-width|auto>
    <associate|project-flag|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?|module_2.tm>>
    <associate|auto-10|<tuple|5|?|module_2.tm>>
    <associate|auto-11|<tuple|1.2.3|?|module_2.tm>>
    <associate|auto-12|<tuple|1.2.4|?|module_2.tm>>
    <associate|auto-13|<tuple|6|?|module_2.tm>>
    <associate|auto-14|<tuple|7|?|module_2.tm>>
    <associate|auto-15|<tuple|1.2.5|?|module_2.tm>>
    <associate|auto-16|<tuple|8|?|module_2.tm>>
    <associate|auto-17|<tuple|1.3|?|module_2.tm>>
    <associate|auto-18|<tuple|1.3.1|?|module_2.tm>>
    <associate|auto-19|<tuple|9|?|module_2.tm>>
    <associate|auto-2|<tuple|1.1|?|module_2.tm>>
    <associate|auto-20|<tuple|1.3.2|?|module_2.tm>>
    <associate|auto-21|<tuple|1.3.3|?|module_2.tm>>
    <associate|auto-22|<tuple|1.3.4|?|module_2.tm>>
    <associate|auto-23|<tuple|1.3.5|?|module_2.tm>>
    <associate|auto-3|<tuple|1|?|module_2.tm>>
    <associate|auto-4|<tuple|1.2|?|module_2.tm>>
    <associate|auto-5|<tuple|1.2.1|?|module_2.tm>>
    <associate|auto-6|<tuple|2|?|module_2.tm>>
    <associate|auto-7|<tuple|3|?|module_2.tm>>
    <associate|auto-8|<tuple|1.2.2|?|module_2.tm>>
    <associate|auto-9|<tuple|4|?|module_2.tm>>
    <associate|fig:bulk_sequencing|<tuple|1|?|module_2.tm>>
    <associate|fig:ccf_density|<tuple|8|?|module_2.tm>>
    <associate|fig:ccf_vs_vaf|<tuple|2|?|module_2.tm>>
    <associate|fig:cn_snv|<tuple|6|?|module_2.tm>>
    <associate|fig:genotype_prior|<tuple|7|?|module_2.tm>>
    <associate|fig:mut_genotype|<tuple|3|?|module_2.tm>>
    <associate|fig:pop_structure|<tuple|4|?|module_2.tm>>
    <associate|fig:pop_structure_simple|<tuple|5|?|module_2.tm>>
    <associate|fig:snv_phylogeny|<tuple|9|?|module_2.tm>>
    <associate|part:module_2.tm|<tuple|?|?|../../../.TeXmacs/texts/scratch/no_name_13.tm>>
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
      mutation.>|<pageref|auto-6>>

      <tuple|normal|<surround|<hidden|<tuple>>||Effect of mutational genotype
      on observed VAF.>|<pageref|auto-7>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of the assumed
      population structure. Here all populations are defined with respect to
      a single mutation. The circular cells are non-malignant and the
      irregularly shaped ones are malignant.>|<pageref|auto-9>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of the
      simplified population structure. In contrast to Figure
      <reference|fig:pop_structure_simple> the mutational genotypes of all
      cells within the reference and variant populations are the same. Note
      the genotypes are different between the populations, as they must be
      since one population has the mutation and the other does
      not.>|<pageref|auto-10>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of allele
      specific copy number profile. The red line is the major copy number and
      the blue line is the minor copy number. The star indicates the location
      of an SNV. The major copy number of this SNV is 2 and the minor is
      1.>|<pageref|auto-13>>

      <tuple|normal|<surround|<hidden|<tuple>>||Schematic of how to illicit
      mutational genotype priors. We assume we have the information from
      Figure <reference|fig:cn_snv> and know the major copy number is 2 and
      the minor copy number is 1. The first two examples correspond to
      mutations which happen prior to the copy number change, hence the total
      copy number of the reference and variant population differ. The third
      example corresponds to the case where the mutation occurs after the
      copy number event. Hence, the copy number of the reference and variant
      population are the same. Furthermore, only a single copy can be mutated
      by the infinite sites assumption.>|<pageref|auto-14>>

      <tuple|normal|<surround|<hidden|<tuple>>||Example posterior densities
      for the cellular prevalence <with|mode|<quote|math>|\<phi\>> computed
      from the model. The top row shows the case for a homozygous diploid
      position (there is a s typo it should be CN=(2,0)). The second row
      shows the posterior for the case illustrated in Figures
      <reference|fig:cn_snv> and <reference|fig:genotype_prior>.>|<pageref|auto-16>>

      <tuple|normal|<surround|<hidden|<tuple>>||Illustration of the
      relationship between evolutionary history and cellular prevalence. On
      the left we have hypothetical evolutionary history, where stars
      indicate mutations and nodes clonal populations. On the write is a
      hypothetical set of cellular prevalence for the
      mutations.>|<pageref|auto-19>>
    </associate>
    <\associate|parts>
      <tuple|module_2.tm|chapter-nr|0|section-nr|0<uninit>|subsection-nr|0>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Inferring
      clonal population structure from SNV data>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Modelling mutation genotype
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|1.2.1<space|2spc>Background
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|1.2.2<space|2spc>Population structure
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|1.2.3<space|2spc>Modelling counts
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|2tab>|1.2.4<space|2spc>Genotype priors
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|1.2.5<space|2spc>The full independent model
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>Clustering mutations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|2tab>|1.3.1<space|2spc>Motivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|2tab>|1.3.2<space|2spc>Mixture models
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|2tab>|1.3.3<space|2spc>Sharing statistical
      strength <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|2tab>|1.3.4<space|2spc>MCMC inference
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|2tab>|1.3.5<space|2spc>Dirichlet process
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>
    </associate>
  </collection>
</auxiliary>