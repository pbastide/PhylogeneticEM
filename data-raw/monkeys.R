## Monkeys dataset
# Taken from: http://www.pnas.org/content/113/8/2158.abstract

## Data : PC1 and PC2 fro brain shape.
# Link: http://www.pnas.org/content/suppl/2016/02/04/1514473113.DCSupplemental/pnas.1514473113.sd01.xls

dat <- read.csv(file = "../data/Aristide2016/pnas.1514473113.sd01.csv")

## Phylogeny
# Link: http://www.pnas.org/content/suppl/2016/02/04/1514473113.DCSupplemental/pnas.1514473113.sd02.pdf

phy <- read.tree(text = "((((((((Saguinus_bicolor:2.998996569,Saguinus_midas:2.998996569)58:6.078720853,Saguinus_mystax:9.077717423)57:1.616541786,Saguinus_fuscicollis:10.69425921)56:4.854715739,((((((Callithrix_jacchus:0.6452483447,Callithrix_penicillata:0.6452483447)64:0.6885388748,Callithrix_kuhlii:1.333787219)63:2.926048241,Callithrix_aurita:4.25983546)62:2.856581709,(((Mico_humeralifer:1.037075975,Mico_chrysoleucus:1.037075975)67:1.037075975,Mico_argentatus:2.07415195)66:3.3464547,Cebuella_pygmaea:5.420606651)65:1.695810519)61:5.82887065,Callimico_goeldii:12.94528782)60:1.573182977,(Leontopithecus_chrysomelas:0.9346547268,Leontopithecus_rosalia:0.9346547268)68:13.58381607)59:1.030504151)55:4.725565576,((Aotus_azarae:1.911998223,Aotus_nigriceps:1.911998223)70:3.549201739,Aotus_vociferans:5.461199962)69:14.81334056)54:0.5657543201,(((Cebus_apella:1.607248907,Cebus_robustus:1.607248907)73:6.847829154,(Cebus_albifrons:3.152380031,Cebus_capucinus:3.152380031)74:5.30269803)72:9.354765376,(Saimiri_sciureus:3.315674187,Saimiri_boliviensis:3.315674187)75:14.49416925)71:3.030451407)53:3.080317337,((((Ateles_geoffroyi:3.232076003,Ateles_marginatus:3.232076003)80:0.641017712,(Ateles_paniscus:2.843114279,Ateles_chamek:2.843114279)81:1.029979435)79:9.185356951,(Lagothrix_lagotricha:9.003321026,(Brachyteles_hypoxanthus:1.766813375,Brachyteles_arachnoides:1.766813375)83:7.236507651)82:4.05512964)78:3.842640595,((((Alouatta_macconnelli:3.07881384,Alouatta_seniculus:3.07881384)87:0.8291311444,Alouatta_caraya:3.907944984)86:1.297523474,(Alouatta_belzebul:4.053668767,Alouatta_guariba:4.053668767)88:1.151799692)85:1.222139909,(Alouatta_palliata:3.651035304,Alouatta_pigra:3.651035304)89:2.776573064)84:10.47348289)77:7.019520919)52:1.399708022,((((Callicebus_moloch:3.799240445,Callicebus_brunneus:3.799240445)93:2.107527241,Callicebus_donacophilus:5.906767685)92:5.76565188,Callicebus_personatus:11.67241957)91:9.045620095,(((Cacajao_calvus:3.405820865,Cacajao_melanocephalus:3.405820865)96:4.09613133,(Chiropotes_satanas:1.876760013,Chiropotes_chiropotes:1.876760013)97:5.625192182)95:5.124611168,((Pithecia_monachus:0.6264591628,Pithecia_irrorata:0.6264591628)99:1.685914079,Pithecia_pithecia:2.312373242)98:10.31419012)94:8.091476297)90:4.602280542)51;")

## Make object
monkeys <- list(phy = phy, dat = dat)