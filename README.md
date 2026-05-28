# transcriptomics jaar 2

# Transcriptomics Analyse: Differentieels Expressieonderzoek naar Reumatoïde Artritis (RA)
## 📁 Inhoudsopgave
1. Introductie & Achtergrond
2. Doel van het Onderzoek
3. Materiaal & Methoden
4. Resultaten
5. Conclusie & Discussie
6. Data Stewardship & Repositorystructuur

## 🧠 Introductie & Achtergrond
Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wereldwijd bij ongeveer 1% van de bevolking voorkomt. De aandoening kenmerkt zich door chronische inflammatie van het gewrichtsslijmvlies (synovium), met name in de handen en voeten. Dit leidt tot progressieve gewrichtsschade, pijn, stijfheid en uiteindelijk tot weefselschade of invaliditeit.
Hoewel de exacte etiologie nog onbekend is, speelt een complexe interactie tussen genetische factoren (zoals specifieke HLA-DRB1-allelen) en omgevingsfactoren (zoals roken) een cruciale rol. Vroege detectie en behandeling zijn essentieel om gewrichtsschade bij tot wel 90% van de patiënten drastisch te beperken. Omdat het ziektebeeld en de symptomen per patiënt sterk variëren, is het identificeren van betrouwbare, vroege biomarkers een grote uitdaging binnen de reumatologie.

## 🎯 Doel van het Onderzoek
Het doel van dit transcriptomics-onderzoek is het identificeren van genen die significant hoger of lager tot expressie komen in het synoviumweefsel van RA-patiënten in vergelijking met gezonde controles. Daarnaast beoogt deze studie inzicht te krijgen in de specifieke biologische processen en metabole/biochemische pathways die door deze expressieverschillen worden ontregeld.

## 🧬 Materiaal & Methoden
### Dataset en Monsters
Voor deze analyse is gebruikgemaakt van RNA-sequencingdata (FASTQ-bestanden) afkomstig uit synoviumbiopten van acht vrouwelijke individuen:
4 x RA-patiënten: Diagnose $> 12$ maanden (established RA), ACPA-positief (anti-CCP).
4 x Gezonde controles: ACPA-negatief.
### Bioinformatische Pipeline (R-pakketten)
De transcriptomics-pipeline is volledig uitgevoerd in R Studio met de volgende core packages:
Data-preparatie & Alignment: Rsubread voor het indexeren van het humane referentiegenoom (hg38/RefSeq GCF_000001405.26) en het uitlijnen van de ruwe reads (align()).
Quantificatie: Rsubread::featureCounts voor het genereren van de count-matrix op exon-niveau met behulp van een bijbehorende GTF-annotatie.
BAM-verwerking: Rsamtools voor het sorteren en indexeren van de gegenereerde alignment-bestanden.
Differentiële Expressie Analyse: DESeq2 voor de statistische bepaling van de genexpressieverschillen tussen de RA- en de controlegroep.
Gen-ID Conversie: org.Hs.eg.db voor de mapping tussen gen-symbolen en ENTREZ ID's.
Functionele Annotatie & Verrijking: * goseq (inclusief geneLenDataBase) voor Gene Ontology (GO) verrijkingsanalyse met ingebouwde correctie voor genlengte-bias.clusterProfiler en KEGGREST voor het identificeren van biologische pathways.
Visualisatie: EnhancedVolcano (volcano plots), pathview (visuele mapping op KEGG-pathwaykaarten) en ggplot2 (GO-barplots).

## 📊 Resultaten

### 1. Differentiële Genexpressie (Volcano Plot)
De analyse met DESeq2 (gevisualiseerd via EnhancedVolcano) toont een duidelijk profiel van genen die significant verschillend gereguleerd zijn tussen RA-patiënten en de controlegroep. Opvallend is dat de volcano-plot in de totale distributie een groter aantal downregulated genen laat zien dan upregulated genen.
De top meest opvallende genen zijn:
Sterkste Overexpressie (Upregulated): MTND5P5 (gekoppeld aan verhoogde inflammatie en RA-risico bij vrouwen), KRT14 (botreparatie/respons op schade), AHSP (bescherming tegen oxidatieve stress), GP9 (bloedplaatjesvorming/megakaryocyten) en IGHV3-53 (productie van auto-antistoffen). Ook SRGN, BCL2A1, COL6A5 en PTGFR vertonen een significant verhoogde expressie.
Sterkste Onderexpressie (Downregulated): HNRNPA3P6 (pseudogen-regulator), CSN2 (beïnvloed door de antioxidant Nrf2), ZNF511-PRAP1 (onderscheid tussen eigen en vreemd DNA), UTP14C (ribosoom-aanmaak) en AFM (vitamine E-transport tegen oxidatieve stress).
### 2. Gene Ontology (GO) Verrijking
De top biologische processen die naar voren komen uit de goseq-analyse zijn direct gerelateerd aan het immuunsysteem. Processen zoals de immuunrespons, inflammatoire reacties en lymfocytenactivatie zijn significant veranderd bij RA-patiënten.
### 3. KEGG Pathway Analyse
Binnen de specifieke reuma-pathway (hsa05323: Rheumatoid arthritis) en gerelateerde signaalroutes (zoals TNF-, IL-17-, en MAPK-signaalroutes) zijn kritieke ontregelingen vastgesteld:
Verhoogde expressie (Upregulated): IL6, IL1β en MMP13. Deze genen drijven de chronische ontsteking aan en zijn direct verantwoordelijk voor kraakbeenschade en gewrichtsvernietiging. Daarnaast is RANK verhoogd, wat leidt tot overactivatie van osteoclasten (botafbraak). Ook chemokinen (CCL-familie) tonen een hogere expressie, wat de migratie van immuuncellen naar het gewricht verklaart.
Opvallende biomarker: Het gen PTGDS (prostaglandine D2 synthase) vertoont een opvallend hoge expressie, wat een sleutelrol suggereert in de chronische ontstekingscascade en mogelijk potentie heeft als vroege biomarker.

## 🛑 Conclusie & Discussie
Dit onderzoek bevestigt dat transcriptomics een krachtige methode is om de moleculaire mechanismen achter Reumatoïde Artritis te ontrafelen. De resultaten tonen een duidelijke verschuiving in het transcriptoom van RA-patiënten, gekenmerkt door de activatie van ontstekingsbevorderende cytokines (IL6, IL1β), botafbrekende factoren (RANK), en het potentiele biomarkargen PTGDS.

### Beperkingen & Aanbevelingen voor Vervolgonderzoek:
1. Steekproefgrootte: De huidige studie is gebaseerd op een kleine populatie ($n=8$). Een grotere steekproef is noodzakelijk om de statistische power te verhogen.
2. Patiëntdiversiteit: De huidige dataset bestaat exclusief uit vrouwelijke patiënten en kent een scheve leeftijdsverdeling tussen de controle- en behandelgroep. Vervolgonderzoek moet ook mannen en leeftijd-gematchte controles bevatten om biologische bias (zoals hormoneffecten) uit te sluiten.
3. Pathway-breedte: De focus lag nu sterk op immuun-gerelateerde (upregulated) pathways. Toekomstig onderzoek moet ook de downregulated pathways nader onderzoeken om te begrijpen welke beschermende mechanismen onderdrukt worden.

## 📁 Data Stewardship & Repositorystructuur

Dit project is ingericht volgens de FAIR-dataprincipes om reproduceerbaarheid en transparantie te garanderen.
Plaintext

├── Assets/                 # Beeldmateriaal, logo's en flowcharts van het werkproces
├── Bronnen/                # Wetenschappelijke literatuur (PDF's) ter ondersteuning van de resultaten
├── Data/
│   ├── Raw/                # Ruwe FASTQ sequencing data (gezipt) en patiëntkenmerken (CSV)
│   └── Processed/          # Gealigneerde data en verwerkte datasets
├── Data_Stewardship/       # Documentatie omtrent databeheer en -opslag conform de richtlijnen
├── Resultaten/             # Gegenereerde plots (Volcano plot, GO-barplots, KEGG pathview-afbeeldingen)
├── Scripts/                # Het volledige, reproduceerbare R-script (Script_R_Casus.R)
└── README.md               # Dit document

Dit project is uitgevoerd als onderdeel van het tweedejaars curriculum Transcriptomics (Jaar 2, Periode 4).