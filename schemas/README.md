# Frictionless Data Schemas

<https://framework.frictionlessdata.io/>

In `igseq.package.yaml` with this directory alongside it:

    name: igseq
    resources:
      - name: subjects
        path: subjects.csv
        schema: schemas/subjects.schema.yaml
      - name: specimens
        path: specimens.csv
        schema: schemas/specimens.schema.yaml
      - name: samples
        path: samples.csv
        schema: schemas/samples.schema.yaml
      - name: runs
        path: runs.csv
        schema: schemas/runs.schema.yaml

And then:

    $ frictionless validate igseq.package.yaml
    ────────────────────────────── Dataset ──────────────────────────────
                        dataset                    
    ┏━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━┓
    ┃ name      ┃ type  ┃ path          ┃ status  ┃
    ┡━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━┩
    │ subjects  │ table │ subjects.csv  │ VALID   │
    │ specimens │ table │ specimens.csv │ INVALID │
    │ samples   │ table │ samples.csv   │ VALID   │
    │ runs      │ table │ runs.csv      │ VALID   │
    └───────────┴───────┴───────────────┴─────────┘
    ...
