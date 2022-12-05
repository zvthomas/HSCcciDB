from django.contrib import admin

from .models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype

admin.site.register(Pathway)
admin.site.register(Receptor)
admin.site.register(Ligand)
admin.site.register(cellClas)
admin.site.register(pathwayAndCelltype)

