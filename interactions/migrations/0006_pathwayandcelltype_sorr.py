# Generated by Django 4.0.2 on 2022-08-17 22:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0005_rename_cellclass_cellclas'),
    ]

    operations = [
        migrations.AddField(
            model_name='pathwayandcelltype',
            name='sorr',
            field=models.CharField(choices=[('s', 'Sending'), ('r', 'Recieving')], default='s', max_length=1),
        ),
    ]