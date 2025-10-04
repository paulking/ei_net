# E-I Net

![E-I Net](https://www.pking.org/research/EINet/EI_Net_logo-600x200.png)

E-I Net is a spiking neural circuit model coded in MATLAB that learns sparse code representations
such as found in the brain. The simulator was developed in 2013.

E-I Net models aspects of the pattern-learning neural networks found the cerebral cortex by using separate
populations of excitatory and inhibitory spiking neurons to find statistical patterns in an input information space.
Using only unsupervised learning (i.e. no training or reference signal), E-I Net learns Gabor-like
receptive fields such as those found in primary visual cortex (area V1).
The spiking circuit simulator is written in MATLAB.

This work was done at the UC Berekeley Redwood Center for Theoretical Neuroscience and
is an evolution of SAILnet by Joel Zylberberg.

**FAQ:** The motivation behind E-I Net is described in the [E-I Net FAQ](https://www.pking.org/research/EINet/ei_net_faq.html).

**Documentation:** The inner workings of E-I Net's spiking circuit simulator is explained in this
[Technical Overview](https://www.pking.org/research/EINet/Neurosim_Technical_Overview.pdf).

**Source code:** The E-I Net source code and the spiking circuit simulator on which it is based (Neurosim)
is written in MATLAB and can be downloaded from GitHub (here) or as a zipfile.

## References

King PD, Zylberberg J, DeWeese MR (2013). Inhibitory interneurons decorrelate excitatory cells to drive sparse code
formation in a spiking model of V1. *J Neuroscience* 33(13): 5475-5485. [abstract](http://www.jneurosci.org/content/33/13/5475.abstract),
[PDF](http://redwood.berkeley.edu/w/images/2/29/King_Zylberberg_DeWeese_E_I_Net_Model_of_V1_JNeurosci_2013.pdf)

King PD, Zylberberg J, DeWeese MR (2012). Inhibitory interneurons enable sparse code formation in a spiking circuit
model of V1. *BMC Neuroscience* 13(Suppl 1):P148. Presented at Computational Neuroscience Society (CNS) 2012.
[PDF](http://www.biomedcentral.com/content/pdf/1471-2202-13-S1-P148.pdf)

Zylberberg J, Murphy JT, DeWeese MR (2011). A sparse coding model with synaptically local plasticity and spiking neurons can account
for the diverse shapes of V1 simple cell receptive fields. *PLoS Computational Biology* 7(10):
e1002250. doi:10.1371/journal.pcbi.1002250.
[article](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002250)

For details of the E-I Net project (2013), see the [EI-Net Home Page](https://www.pking.org/research/EINet/index.html)

