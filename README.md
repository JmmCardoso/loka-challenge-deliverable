# GeneXOmics Cloud Migration & Automation

## Project Overview

This repository contains the final deliverable for the Loka Technical Challenge.
It presents the technical implementation and architectural design for migrating GeneXOmics’ Perturb-Seq data workflows from a manual, on-premise/hybrid setup to a fully automated, scalable, and event-driven cloud infrastructure on Amazon Web Services (AWS).

---

## Repository Structure & Deliverables for Loka Challenge

### [1. Bioinformatics Pipeline](./pipeline/)

* **Location**: `/pipeline`
* **Description**: A Nextflow (DSL2) pipeline encapsulating `bcl2fastq` and `cellranger`.
* **Features**:
  * Dockerized execution modules.
  * Automatic cell count detection logic.
  * Setup scripts for local testing.

### [2. Platform Analysis](Platform_Analysis.md)

* **Document**: `Platform_Analysis.md`
* **Context**: Technical evaluation of cloud platforms.
* **Recommendation**: **AWS Batch** (Managed) due to Spot Instance savings (>50%) and native Nextflow integration.

### [3. Architecture Diagram](Architecture_Diagram_Support.md)

* **Diagram**: `Step3_AWS_Diagram.png`
* **Support Document**: `Architecture_Diagram_Support.md`

### [4. Cost Estimate](AWS_Cost_analysis.pdf)

* **Document**: `AWS_Cost_Analysis.pdf`
* **Context**: Monthly operational cost projection.
* **Highlight**: Estimated run cost of **<$200/month** for infrastructure, leveraging Spot pricing.

---

*For GeneXOmics Infrastructure Review*
