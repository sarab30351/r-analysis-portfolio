# Some background information that supports the rationale of my study:

Breast cancer can be split into molecular "subtypes" (like those caught by the PAM50 test + Claudin-low, which PAM50 doesn't usually capture) that help predict how patients will do. But how much extra risk each subtype carries isn't always constant, it has been shown to change over the years, and for some patients the subtype may matter less and less the longer they survive.

Most past studies accounted for this by adjusting for standard clinical and pathology factors, or by using existing risk-scoring tools. And the long-term studies that used the METABRIC dataset mostly looked at a different kind of genomic classification, or only looked at specific smaller groups of patients.

What's still unclear (ad why I chose to investigate it) is this: if you look at PAM50/Claudin-low subtypes on top of an already established risk tool (the Nottingham Prognostic Index) and age, and you allow for the fact that each subtype's effect on risk can rise or fall at different points during follow-up — does subtype still add meaningful prognostic information? That question hasn't been fully answered yet (atleast to my knowledge and research).
I also wanted to give a prediction element to this project, but I'm still contemplating about it (mostly because of the limited applicability; it would be great for practice but it would barely add any value).

Talking about research, here's some literature on this topic, in chronological order of developments:

1. Parker and colleagues establish PAM50's independent prognostic value -> https://pubmed.ncbi.nlm.nih.gov/19204204/
2. Nielsen and colleagues demonstrate that the subtype effect differs between early and later survival -> https://pmc.ncbi.nlm.nih.gov/articles/PMC2970720/
3. Caan and colleagues demonstrate the phenomenon explicitly for recurrence and mortality and show just how strongly the hazard ratios change over time -> https://aacrjournals.org/cebp/article/23/5/725/70359/
4. Rueda and colleagues establish long-term molecular recurrence dynamics specifically in METABRIC cohort -> https://www.nature.com/articles/s41586-019-1007-8
5. Lundgren and colleagues provide a modern long follow-up PAM50 analysis in which non-proportionality is explicitly recognized -> https://www.nature.com/articles/s41523-022-00423-z
6. Richman and colleagues show, using METABRIC itself, that the prognostic importance of intrinsic subtype may largely diminish by very late follow-up -> https://pubmed.ncbi.nlm.nih.gov/38709373/
7. Zhen and colleagues demonstrated the strong correlation between NPI and molecular subtype, leading me to use NPI as the comparator -> https://pubmed.ncbi.nlm.nih.gov/29088770/
