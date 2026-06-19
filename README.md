 <img width="1699" height="926" alt="image" src="https://github.com/user-attachments/assets/9fd35771-4cc1-4ebd-9f39-6909a61d93e2" />


📁 Inhoud
---
★  `data/`

★ `scripts/`

★ `resultaten/`

★ `referenties/`

📖 Inleiding 
---
Reumatoïde Artritis (RA) is een auto-immuunziekte die gekenmerkt wordt door ontstekingen aan het gewrichtsslijmvlies, deze ontstekingen ontstaan doordat het immuunsysteem de eigen gewrichten aanvalt. Hoewel de precieze oorzaak nog onbekend is, is het wel duidelijk dat er zowel genetische als omgevingsfactoren een rol spelen. Zoals een polymorfisme, een veelvoorkomende mutatie in DNA bij verschillende mensen, in het HLA-DRB1-regio van het major histocompatibility complex (MHC). HLA-DRB1-gen speelt een rol bij het presenteren van antigenen aan T-cellen, sommige varianten van dit gen vergroten de kans dat het immuunsysteem lichaamseigen eiwitten als vreemd herkennen (Padyukov L. 2022). Daarnaast spelen ook omgevingsfactoren zoals luchtwegontstekingen en ernstige tandvleesontsteking (chronische parodontitis) veroorzaakt door *Porphyromonas gingivalis* een rol bij het ontstaan en progressie van RA. Antistoffen en ontstekingsmarkers zijn voordat er echte klachten zijn al detecteerbaar, de bekenste serologische markers voor RA zijn de reumafactor (RF) en de anti-citrullinated peptide/proteins antibodies (ACPA) (Kolarz et al., 2018). Door deze ontregelde immuunrespons waarbij T- en B-lymfocyten betrokken zijn, kunnen veranderingen in genexpressie bijdragen aan een beter inzicht in de biologische processen die betrokken zijn bij het ontstaan en de progressie van RA. Een transcriptomics analyse met R wordt gebruikt om de betrokken genen en pathways te identificeren. Het doel van dit onderzoek is om een beter inzicht te krijgen welke genen en pathways zijn betrokken bij RA.

🧬 Methode
---
<p align="center">
  <img width="750" height="350" alt="image" src="https://github.com/user-attachments/assets/2b2f14da-a67f-4075-a643-077ff2d002c3" />
</p>

 <p align="center">
 <em><strong>Figuur 1.</strong> Workflow voor de Transcriptomics analyse van de genen en pathways bij RA patiënten en controle groep.</em>


Om te onderzoeken welke genen en pathways zijn betrokken bij RA. Is er een biopt genomen uit het gewrichtsslijmvlies van 4 vrouwen met RA en 4 vrouwen zonder RA (controle groep). De vrouwen met RA testen positief op ACPA en de vrouwen zonder negatief, de data zijn afkomstig uit het onderzoek van Platzer et al. (2019). Allereerst werd er met R (versie 4.5.2) het humane referentiegenoom GRCh38.p14 (accesion number: GCF_000001405.40) geïndexeerd met behulp van het Rsubread package (versie 2.24.0), om het alignen snel en gemakkelijk te laten verlopen. Met behulp van de package Rsubread (versie 2.24.0) is hierna een count matrix gemaakt, om te achterhalen hoeveel reads er per gen zijn gemapt. Aan de hand van deze matrix is de differentiële expressie-analyse uitgevoerd met de package DESeq2 (versie 1.50.2). Tijdens de analyse werd de genexpressie tussen de vrouwen met RA en de controle groep met elkaar vergeleken, om te identificeren welke genen een significant verhoogde of verlaagde expressie toonden. Om dit verschil te visualiseren is er een Volcano plot gemaakt met behulp van de package EnhancedVolcano (versie 1.28.2). Voor het in kaart brengen van de biologische processen die betrokken zijn bij de genen die significant een veranderde expressie toonde is er een Gene Ontology (GO)-analyse uitgevoerd. Hiervoor is gebruik gemaakt van de packages clusterProfiler (versie 4.18.4 ). Als laatste is de KEGG-pathway analyse gebruikt om de significant meer aanwezige pathways te identificeren met behulp van de package Pathview (versie 1.50.0)
 
📊 Resultaten 
---
#### *Volcano Plot*
Na het uitvoeren van de differentiële expressie-analyse tussen de RA- en controlesamples op basis van de count matrix, is hier een Volcano plot van gemaakt [Figuur 2](resultaten/Figuur%201.%20Volcano%20plot.pdf). Hierin is de Log<sub>2</sub>FoldChange uitgezet tegen de -Log<sub>10</sub>*P*, de verticale stippellijnen geeft de drempel voor de log<sub>2</sub>-fold change (±1) weer en de horizontale stippellijn de significantiedrempel (p = 0,05). In totaal zijn er 29.407 genen geanalyseerd, hiervan zijn er 4.572 genen die een significant verschil toonde tussen de RA- en controlesamples. Deze 4.572 genen worden verder gebruikt voor de GO analyse en KEGG pathway analyse. 

<p align="center">
  <img width="500" alt="Volcano plot" src="https://github.com/user-attachments/assets/cfadfa68-9207-4063-838b-7db719d197ca">
</p>

<p align="center">
  <em><strong>Figuur 2.</strong> De x-as geeft de log<sub>2</sub>-fold change weer en de y-as de −log<sub>10</sub>(p-waarde). Genen met een significante aangepaste p-waarde (padj &lt; 0,05) en |log<sub>2</sub>FC| &gt; 1 zijn weergegeven in rood.</em>
</p>


---
#### *GO Analyse* 
In [Figuur 3](resultaten/Figuur%202.%20DotPlot%20GO%20analyse.png) is de Dotplot van de Gene Ontology (GO) analyse van de 4.572 genen die significant differentiële genexpressie toonden. Op de x-as staat de verhouding van het aantal genen uit de dataset dat geassocieerd is met de specifieke GO-term ten opzichte van het totale aantal genen binnen de analyse(GeneRatio). De grootte van de punten geeft het aantal genen (Count) weer en de kleur geeft de aangepaste p-waarde (p.adjust) aan. De analyse laat zien dat met name processen die te maken hebben met de immuunrespons significant verrijkt zijn, zoals de lymphocyte differentiation, leukocyte mediated immunity, T cell differentiation en B cell activation.

<div align="center">
  <img src="https://github.com/user-attachments/assets/6237665b-1d11-4a45-b34c-40599292d7e8" width="500"/>

  *<p><b>Figuur 3. Dotplot van de GO-analyse.</b> De 15 significant verrijkte biologische processen zijn weergegeven. De grootte van de bollen geeft het aantal genen (Count) weer, de kleur de aangepaste p-waarde (p.adjust) en de x-as de GeneRatio.</p>*
</div>

---
#### *KEGG-pathway analyse*
Om inzicht te krijgen in de de biologische pathways waarin de genen die significant differentiële genexpressie toonden betrokken zijn, werd er een KEGG-pathway analyse uitgevoerd. In [Figuur 4](resultaten/Figuur%202.%20KEGG%Dotplot%.png)

<div style="display: flex; gap: 10px;">
  <img src="https://github.com/user-attachments/assets/0da484f9-9a2d-45bb-bc1a-52167f586a13" width="49%" />
  <img src="https://github.com/user-attachments/assets/f21e78fc-e298-46e5-aaab-e62368ffa819" width="49%" />
</div>

---
#### *KEGG-pathway*
<div style="display: flex; gap: 10px;">
  <img src="https://github.com/user-attachments/assets/d1b6df89-a5b5-4da7-bf5a-b518287adac4" width="49%" />
  <img src="https://github.com/user-attachments/assets/bcc41cb4-35b1-4295-b7ed-fc326aebe966" width="49%" />
</div>


📌 Conclusie
---
