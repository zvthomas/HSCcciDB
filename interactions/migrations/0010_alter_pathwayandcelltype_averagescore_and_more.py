# Generated by Django 4.0.2 on 2022-08-18 22:46

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0009_pathwayandcelltype_averagescore'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pathwayandcelltype',
            name='averageScore',
            field=models.FloatField(default=0),
        ),
        migrations.AlterField(
            model_name='pathwayandcelltype',
            name='hscPercent',
            field=models.FloatField(default=0),
        ),
    ]
