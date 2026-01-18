# Platform Execution Analysis

## Cloud Execution Platform: AWS Batch

> **Methodology Note**: This analysis is based on a review of AWS documentation, Nextflow executor specifications, and published genomics workflow best practices. References are provided throughout to support each recommendation.

---

## 1. Platform Selection

Three cloud execution platforms were evaluated for deploying the Perturb-Seq pipeline: **AWS Batch**, **AWS HealthOmics**, and **Terra.bio**.

| Criteria | AWS Batch | AWS HealthOmics | Terra.bio |
|----------|-----------|-----------------|-----------|
| **Workflow Language** | Nextflow (native) | WDL/CWL only | WDL primary |
| **Migration Effort** | None | High (rewrite) | High (rewrite) |
| **Cost Model** | Pay-per-use + Spot | Managed pricing | Cloud + platform fee |
| **Infrastructure Control** | Full | Limited | Limited |
| **AWS Integration** | Native | Native | Multi-cloud (GCP-first) |

**Recommendation**: **AWS Batch** is the recommended platform for GeneXOmics.

**Rationale**:
1. **Native Nextflow support** — AWS Batch is a first-class executor in Nextflow, meaning the pipeline developed in Step 1 requires no code modifications, only configuration changes [1].
2. **Cost optimization** — Supports EC2 Spot Instances, which AWS documentation indicates can reduce compute costs by 60–90% compared to On-Demand pricing [2].
3. **AWS ecosystem alignment** — GeneXOmics is migrating from GCP to AWS; Batch integrates natively with S3, Lambda, and CloudWatch, and is compatible with Quilt's S3-backed data management [3].

---

## 2. Setup Requirements

Based on AWS Batch documentation and Nextflow's executor specifications [1][4], the following infrastructure components are required:

### 2.1 Compute Environment

- **Type**: Managed Compute Environment (AWS handles EC2 provisioning)
- **Instance Strategy**: Mixed instance types with Spot capacity
  - `c5` family (compute-optimized) for alignment tasks
  - `r5` family (memory-optimized) for large sample processing
- **Scaling**: Auto-scaling from 0 to maximum vCPUs based on queue depth

### 2.2 Storage Layer (Amazon S3)

S3 serves as the shared data layer for the pipeline [5]:

| Bucket Purpose | Lifecycle Policy |
|----------------|------------------|
| `raw-data/` | Transition to Glacier after 30 days |
| `work/` (Nextflow scratch) | Delete after 7 days |
| `results/` | Versioning enabled, long-term retention |
| `references/` | Intelligent-Tiering |

### 2.3 Container Registry (Amazon ECR)

Docker images for pipeline tools (`bcl2fastq`, `cellranger`, `scanpy`) are stored in ECR. AWS recommends ECR over Docker Hub for production workloads due to rate-limiting, security, and latency benefits [6].

### 2.4 IAM & Networking

- **IAM Roles**: Separate roles for Batch service, EC2 instances, and Nextflow submission (least-privilege principle) [7]
- **VPC Configuration**: Compute instances run in private subnets; VPC Endpoints for S3 and ECR eliminate NAT Gateway costs and keep traffic internal to AWS [8]

### 2.5 Nextflow Configuration

The key configuration parameters for AWS Batch execution [1]:

```groovy
process.executor = 'awsbatch'
process.queue = 'genexomics-queue'
workDir = 's3://genexomics-genomics/work'
aws.region = 'us-east-1'
```

No changes to pipeline logic are required—only configuration.

---

## 3. Execution Approach

### 3.1 Job Submission

Nextflow translates each pipeline process into an AWS Batch job. AWS Batch dynamically provisions EC2 instances, pulls containers from ECR, and executes tasks in parallel [4].

### 3.2 Automation (Recommended for Production)

To meet GeneXOmics's goal of reduced manual intervention, an event-driven trigger is recommended:

1. **AWS DataSync** transfers data from on-premise NAS to S3 [9]
2. **S3 Event Notification** detects the presence of `RTAComplete.txt`
3. **EventBridge → Lambda** submits the AWS Batch job automatically [10]

This pattern eliminates the need for manual pipeline submission once sequencing completes.

### 3.3 Monitoring

| Tool | Purpose |
|------|---------|
| AWS Batch Console | Job status, retries, exit codes |
| CloudWatch Logs | Container stdout/stderr for debugging |
| CloudWatch Alarms | Alerts on job failures or queue backlog |
| Seqera Platform (optional) | Web UI for pipeline visualization [11] |

---

## 4. Advantages for GeneXOmics

| Advantage | Relevance to GeneXOmics |
|-----------|-------------------------|
| **Native Nextflow support** | No pipeline rewrite; existing code works directly |
| **Spot Instance pricing** | 60–90% compute cost reduction for bursty workloads (8 runs/week) |
| **Scale-to-zero** | No cost when idle; ideal for SMB budget constraints |
| **Parallel execution** | Multiple samples processed simultaneously |
| **Automatic retries** | Self-healing on transient failures or Spot interruptions |
| **AWS integration** | Clean fit with S3, Lambda, Quilt, and planned infrastructure |

---

## 5. Limitations and Mitigations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| **Cold start latency** | 2–5 min to provision first EC2 instance | Acceptable for batch genomics; optional warm pool if needed [4] |
| **Spot interruptions** | Instances can be reclaimed by AWS | Nextflow checkpointing enables resume; On-Demand fallback available [12] |
| **Operational complexity** | More components than fully managed platforms | Infrastructure-as-Code (Terraform/CloudFormation) simplifies deployment [13] |
| **No genomics-specific abstractions** | Unlike HealthOmics, Batch is generic | Flexibility is advantageous for custom Perturb-Seq workflows |

---

## 6. Performance Considerations

### 6.1 Expected Execution Times

Based on Cell Ranger documentation and typical Perturb-Seq workloads [14]:

| Pipeline Stage | Estimated Runtime (Parallel) |
|----------------|------------------------------|
| BCL → FASTQ | ~10 min |
| FastQC | ~2 min |
| Cell Ranger count | ~8–12 min |
| Downstream analysis | ~3–5 min |
| **Total (per run)** | **~20–30 min** |

*Note: Sequential processing of the same workload would take several hours.*

### 6.2 Scalability

- **Horizontal**: AWS Batch scales linearly with sample count; no architectural changes needed as throughput increases [4]
- **Vertical**: Larger instance types (e.g., `r5.4xlarge`) available for high-cell-count samples
- **Future-proof**: Suitable for reprocessing historical datasets or sudden increases in sequencing volume

---

## 7. Comparison with Alternative Platforms

### AWS HealthOmics
- **Pros**: Purpose-built for genomics; built-in HIPAA compliance
- **Cons**: Requires WDL—complete pipeline rewrite; managed pricing adds overhead
- **Verdict**: Better suited for clinical/regulated workflows, not GeneXOmics's research focus [15]

### Terra.bio
- **Pros**: Strong academic community; collaborative workspace features
- **Cons**: GCP-first architecture (misaligned with AWS migration); WDL-focused
- **Verdict**: Better suited for multi-institutional research collaborations [16]

---

## 8. Summary

AWS Batch is the recommended execution platform for GeneXOmics because it:

- Executes the existing Nextflow pipeline without modification
- Reduces compute costs by 60–90% through Spot Instances
- Scales automatically with workload demand
- Integrates natively with the planned AWS infrastructure
- Supports full automation from sequencer to results

This recommendation aligns with GeneXOmics's stated goals of increased reliability, reduced human intervention, and fast processing times.

---

## References

1. Nextflow Documentation — AWS Batch Executor: https://www.nextflow.io/docs/latest/aws.html
2. AWS EC2 Spot Instances — Pricing: https://aws.amazon.com/ec2/spot/pricing/
3. AWS Batch Overview: https://aws.amazon.com/batch/
4. AWS Batch User Guide — Compute Environments: https://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html
5. Amazon S3 — Lifecycle Configuration: https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lifecycle-mgmt.html
6. Amazon ECR — Best Practices: https://docs.aws.amazon.com/AmazonECR/latest/userguide/best-practices.html
7. AWS IAM — Security Best Practices: https://docs.aws.amazon.com/IAM/latest/UserGuide/best-practices.html
8. AWS VPC Endpoints: https://docs.aws.amazon.com/vpc/latest/privatelink/vpc-endpoints.html
9. AWS DataSync: https://aws.amazon.com/datasync/
10. Amazon EventBridge — Event-Driven Architecture: https://docs.aws.amazon.com/eventbridge/latest/userguide/eb-what-is.html
11. Seqera Platform (Nextflow Tower): https://seqera.io/platform/
12. AWS Batch — Spot Instance Best Practices: https://docs.aws.amazon.com/batch/latest/userguide/spot_fleet_bes_practice.html
13. AWS CloudFormation: https://aws.amazon.com/cloudformation/
14. 10x Genomics — Cell Ranger System Requirements: https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-system-requirements
15. AWS HealthOmics: https://aws.amazon.com/healthomics/
16. Terra.bio Platform: https://terra.bio/
