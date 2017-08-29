mybiom
prune_taxa(rowSums(otu_table(mybiom))>5,mybiom)
prune_taxa(rowSums(otu_table(prune_samples(sample_data(mybiom)[[11]]=="Dessert",mybiom)))>0,prune_samples(sample_data(mybiom)[[11]]=="Dessert",mybiom))
prune_taxa(rowSums(otu_table(prune_samples(sample_data(mybiom)[[11]]=="Cider",mybiom)))>0,prune_samples(sample_data(mybiom)[[11]]=="Cider",mybiom))
prune_taxa(rowSums(otu_table(prune_samples(sample_data(mybiom)[[11]]=="Dessert",mybiom)))>5,prune_samples(sample_data(mybiom)[[11]]=="Dessert",mybiom))
prune_taxa(rowSums(otu_table(prune_samples(sample_data(mybiom)[[11]]=="Cider",mybiom)))>5,prune_samples(sample_data(mybiom)[[11]]=="Cider",mybiom))

sum(counts(dds,normalized=T))/286
sum(counts(dds[dds$orchard=="Dessert"],normalized=T))/144
sum(counts(dds[dds$orchard=="Cider"],normalized=T))/142
