
The Challenge: 

## **Introduction**

GeneXOmics, a rapidly growing genomics company, has reached out to Loka for assistance in addressing critical inefficiencies in their data management and cloud infrastructure. Currently, their highly manual workflow results in significant delays and requires considerable human intervention. GeneXOmics’s account manager from AWS provided Loka with some notes. Read them carefully.

**Company Overview**

GeneXOmics is a rapidly growing genomics company in the SMB segment, having doubled in size over the past year. The company primarily operates on Google Cloud Platform (GCP) but has recently started using Quilt for data management, which has led to a need for Amazon S3 integration. GeneXOmics aims to migrate from GCP to AWS and build a scalable and future-proof infrastructure on AWS.

Their immediate priority is to migrate their Perturb-Seq data to the cloud in a scalable manner, addressing current inefficiencies. Their secondary goal is to integrate metadata from Benchling and SmartSheets into the Quilt Packages generated for each run.

**GeneXOmics current workflow**:

- Sequencing Process: The current workflow involves sequencing data being transferred from 10x Genomics sequencers to a NAS and then manually copied into Quilt. This process is highly manual and inefficient, with only 50% of sequencing runs completed without the need for manual intervention. Two full-time employees dedicate half of their time to ensuring the process runs smoothly, suggesting potential architectural issues on GCP.
- Equipment and Data: The company operates four sequencers, which handle approximately eight sequencing runs per week, generating around 140GB to 190GB of data weekly, equating to about 1.1TB of data per month.

**Key Goals for Sequencing Migration:**

1. Increased Robustness and Reliability: Ensuring that the sequencing process is more reliable, reducing the need for manual intervention.
2. Reduced Human Involvement: Automating processes to minimize the number of people needed to manage sequencing runs.
3. Fast Processing Times: Improving the infrastructure to enable fast parallel processing of sequencing data.

## **The Challenge**

Use the information above to complete the challenge to the best of your ability. There are no right or wrong answers, we’re looking to understand your thought process about this client in particular and have a discussion with you about it afterwards. In this technical assessment, we place great value on the candidate's research skills, recognizing that the ability to effectively find and apply relevant information is crucial in solving complex problems.

**1. Simple pipelining: Code and Documentation**

Using freely accessible raw data in standard formats, implement a simple bioinformatics pipeline using a workflow orchestrator like **Nextflow** or **WDL**, focusing on basic genomic data processing. Develop a simple pipeline that includes the Primary and Secondary stages of analysis for the same type of data that GeneXOmics has.

Create a Git repository to upload the pipeline code. The repository must include a README.MD file containing the following:

- Instructions for installing dependencies.
- Example inputs and outputs of the pipeline.
- Step-by-step guide on how to run the pipeline, with clear examples.

**Documentation**:

- The documentation should be clear and concise, explaining how to install the required tools and execute the pipeline.
- Include a small test dataset to validate the pipeline.

**2.  Platform Execution Analysis**

Choose ONE cloud execution platform where you would deploy and run your pipeline (eg. AWS HealthOmics, Terra, AWS Batch…)

For the platform you choose, provide:

- **Setup requirements**: What configuration is needed to run your pipeline
- **Execution approach**: How you would submit and monitor runs
- **Advantages**: Key benefits of this platform for GeneXOmics's use case
- **Limitations**: Constraints or challenges specific to this platform
- **Performance considerations**: Expected execution time and scalability

**Note**: Hands-on experience with these platforms is not required. If you have practical experience with the platform you choose, please share specific implementation details. Otherwise, provide a well-researched analysis based on documentation and best practices

**3. Architecture Diagram + Architecture component description**

Prepare an architecture diagram with the AWS Services needed to implement an automated genomics pipeline infrastructure for GeneXOmics, as well as the ones used to keep the data organized and accessible. The architecture diagram should be supported by a legend that details the architecture components and steps. Use Miro, Draw.io or similar for this segment.

**4. Cost estimates**

Using the AWS Pricing Calculator [**https://calculator.aws/#/estimate**](https://calculator.aws/#/estimate), provide a cost estimate of running the workflows using the pipelines chosen above and storing the data in S3 or other services of your choosing. You don’t need to add the cost of running Quilt in your calculations.
