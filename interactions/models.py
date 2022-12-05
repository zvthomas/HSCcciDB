from django.db import models
from django.urls import reverse

# Create your models here.


class Pathway(models.Model):
    """Model representing a biological pathway"""
    
    PATH_STATUS = (
        ('c', 'Cell-Cell Contact'),
        ('e', 'ECM-Receptor'),
        ('s', 'Secreted Signaling')
    )

    name = models.CharField(max_length = 10)
    interaction_type = models.CharField(max_length = 1, choices = PATH_STATUS, default = 's')
    hscPercent = models.FloatField()

    ligands = models.ManyToManyField('Ligand', related_name='pathway_ligand')
    receptors = models.ManyToManyField('Receptor', related_name='pathway_receptor')

    evidences = models.TextField(default = '')

    class Meta:
        ordering = ['interaction_type', 'name']

    def __str__(self):
        """String for representing a Model object."""
        return self.name
        
    @property
    def get_absolute_url(self):
        return reverse('pathway-detail', args=[str(self.id)])

    def get_pathway_number(self):
        co = pathwayAndCelltype.objects.filter(pathway = self).filter(hscPercent__gt=0.05).count()
        return int(co)
    
    def get_highest_interaction_score(self):
        co = pathwayAndCelltype.objects.filter(pathway = self).filter(hscPercent__gt=0.05)
        maxi = 0
        for c in co:
            if c.averageScore > maxi:
                maxi = c.averageScore
        return float(maxi)

class cellClas(models.Model):
    
    CELL_STATUS = (
        ('n', 'Non-hematopoietic'),
        ('i', 'Blood & Immune Cells'),
        ('p', 'HPC')
    )
    
    name = models.CharField(max_length=50)
    sendingPathways = models.ManyToManyField(Pathway, related_name='sendingPathways')
    receivingPathways = models.ManyToManyField(Pathway, related_name='receivingPathways')

    cell_type = models.CharField(max_length = 1, choices = CELL_STATUS, default = 'n')

    class Meta:
        ordering = ['cell_type', 'name']

    def __str__(self):
        return self.name

    @property
    def get_absolute_url(self):
        return reverse('cellclass-detail', args=[str(self.id)])

    def get_pathway_number(self):
        co = pathwayAndCelltype.objects.filter(celltype = self).filter(hscPercent__gt=0.05).count()
        return int(co)
    
    def get_highest_interaction_score(self):
        co = pathwayAndCelltype.objects.filter(celltype = self).filter(hscPercent__gt=0.05)
        maxi = 0
        for c in co:
            if c.averageScore > maxi:
                maxi = c.averageScore
        return float(maxi)

class Ligand(models.Model):
    name = models.CharField(max_length = 10)
    receptors = models.ManyToManyField('Receptor')
    pathways = models.ManyToManyField(Pathway)
    evidences = models.TextField(default = '')

    class Meta:
        ordering = ['name']
    
    def __str__(self):
        return self.name
    
    def get_absolute_url(self):
        return reverse('ligand-detail', args=[str(self.id)])

class Receptor(models.Model):
    name = models.CharField(max_length = 10)
    ligands = models.ManyToManyField(Ligand)
    pathways = models.ManyToManyField(Pathway)
    evidences = models.TextField(default = '')
    
    class Meta:
        ordering = ['name']

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('receptor-detail', args=[str(self.id)])

class pathwayAndCelltype(models.Model):
    
    SIGNAL_STATUS = (
        ('s', 'Sending'),
        ('r', 'Receiving')
    )

    pathway = models.ForeignKey(Pathway,on_delete=models.SET_NULL, null=True)
    celltype = models.ForeignKey(cellClas,on_delete=models.SET_NULL, null=True)
    hscPercent = models.FloatField(default = 0)
    averageScore = models.FloatField(default = 0)
    
    sorr = models.CharField(max_length = 1, choices = SIGNAL_STATUS, default = 's')
    
    class Meta:
        ordering = ['sorr', 'pathway', 'celltype', '-hscPercent', '-averageScore']

    def __str__(self):
        return self.pathway.name + "_" + self.celltype.name + "_" + self.sorr

class pathwayCorrelations(models.Model):
    
    SIGNAL_STATUS = (
        ('s', 'Sending'),
        ('r', 'Receiving')
    )

    pathway1 = models.ForeignKey(Pathway,on_delete=models.SET_NULL, null=True, related_name='pathway1')
    p1a = models.CharField(max_length = 1, choices = SIGNAL_STATUS, default = 's')
    pathway2 = models.ForeignKey(Pathway,on_delete=models.SET_NULL, null=True, related_name='pathway2')
    p2a = models.CharField(max_length = 1, choices = SIGNAL_STATUS, default = 's')
    correlation = models.FloatField(default = 0)
    pval = models.FloatField(default = 0)