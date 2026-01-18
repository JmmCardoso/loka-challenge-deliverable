# Architecture Diagram: GeneXOmics Automated Pipeline Infrastructure

> **Note**: This document provides the component legend and data flow description to accompany the architecture diagram (see `Step3_AWS_Diagram.png` for the visual diagram).

---

## Overview

This architecture enables GeneXOmics to transition from a manual, on-premise workflow to a fully automated, event-driven pipeline on AWS. The design executes the Nextflow pipeline developed in Step 1 without modification.

**Key Design Principles**:

- **Event-driven**: Pipeline triggers automatically when sequencing completes
- **Scale-to-zero**: No compute cost when idle; resources provisioned on demand
- **Integrated**: Metadata from Benchling/SmartSheets flows into Quilt packages automatically

---

## Data Flow Steps

### Step 1: Automated Ingestion

| Component | Service | Role |
|-----------|---------|------|
| Sequencer | On-premise | Writes BCL data to Local NAS |
| DataSync Agent | On-premise VM | Monitors NAS, initiates transfer |
| AWS DataSync | AWS | Securely transfers data to S3 via TLS |
| S3 Raw Data | AWS | Landing zone for incoming sequencing runs |

**Flow**: Sequencer → NAS → DataSync Agent → AWS DataSync → S3 Raw Bucket

### Step 2: Event-Driven Trigger

| Component | Service | Role |
|-----------|---------|------|
| S3 Event | AWS | Emits `ObjectCreated` when `RTAComplete.txt` uploads |
| EventBridge | AWS | Matches event rule, routes to Lambda |
| Lambda (Trigger) | AWS | Parses event, submits job to AWS Batch |
| DynamoDB | AWS | Logs run status ("Started", "Processing", "Complete") |

**Flow**: S3 Event → EventBridge Rule → Lambda → AWS Batch + DynamoDB

### Step 3: Pipeline Execution

| Component | Service | Role |
|-----------|---------|------|
| AWS Batch | AWS | Job scheduler; provisions compute dynamically |
| Job Queue | AWS Batch | Queues pipeline jobs by priority |
| Compute Environment | AWS Batch | Manages EC2 Spot fleet (c5/r5 instances) |
| ECR | AWS | Stores Docker images (cellranger, bcl2fastq, scanpy) |
| Nextflow Head Node | EC2 | Orchestrates pipeline; submits child jobs |
| Worker Nodes | EC2 Spot | Execute individual pipeline tasks |

**Flow**: Batch → Provisions EC2 → Nextflow orchestrates → Workers execute tasks

### Step 4: Storage Layer

| Bucket | Purpose | Lifecycle Policy |
|--------|---------|------------------|
| `s3://genexomics-raw` | Immutable raw BCL/FASTQ | Glacier after 30 days |
| `s3://genexomics-work` | Nextflow scratch directory | Delete after 7 days |
| `s3://genexomics-references` | Genome indices (GRCh38) | Intelligent-Tiering |
| `s3://genexomics-results` | Final outputs (matrices, reports) | Versioned, long-term |

### Step 5: Post-Processing & Integration

| Component | Service | Role |
|-----------|---------|------|
| Lambda (Packager) | AWS | Triggers on pipeline completion |
| Benchling API | External | Source of experimental metadata (LIMS) |
| SmartSheets API | External | Source of project tracking data |
| Quilt Registry | External (S3-backed) | Immutable, searchable data packages |

**Flow**: S3 Results → Lambda → Fetch metadata → Create Quilt Package

### Step 6: Monitoring & Alerting

| Component | Service | Role |
|-----------|---------|------|
| CloudWatch Logs | AWS | Container stdout/stderr for debugging |
| CloudWatch Alarms | AWS | Triggers on job failures or anomalies |
| SNS | AWS | Publishes alerts to Slack/Email |

---

## Component Legend

### Compute Services

| Service | Description |
|---------|-------------|
| **AWS Batch** | Managed job scheduler; dynamically provisions EC2 based on queue depth |
| **EC2 Spot Instances** | Worker nodes at 60-90% cost reduction vs On-Demand |
| **AWS Lambda** | Serverless functions for triggering pipelines and post-processing |
| **Amazon ECR** | Private Docker registry for pipeline containers |

### Storage Services

| Service | Description |
|---------|-------------|
| **Amazon S3** | Object storage for all data (raw, intermediate, results) |
| **Amazon DynamoDB** | NoSQL database for tracking run status |

### Integration Services

| Service | Description |
|---------|-------------|
| **AWS DataSync** | Automated, encrypted data transfer from on-premise NAS |
| **Amazon EventBridge** | Event router; triggers Lambda on S3 events |
| **Amazon SNS** | Notification service for alerts (Slack, Email) |
| **Amazon CloudWatch** | Centralized logging and monitoring |

### Networking & Security

| Service | Description |
|---------|-------------|
| **Amazon VPC** | Private network; compute runs in isolated subnets |
| **VPC Endpoints** | Private S3/ECR access without internet traversal |
| **AWS KMS** | Encryption keys for data at rest |
| **IAM Roles** | Least-privilege access control per service |

---

## Security Features

1. **Encryption at Rest**: S3 buckets, DynamoDB, EBS volumes encrypted via AWS KMS
2. **Encryption in Transit**: All transfers use TLS 1.2+
3. **Network Isolation**: Compute in private subnets; no public internet exposure
4. **Least Privilege**: IAM roles scoped to specific actions (e.g., Batch role cannot delete logs)

---

## References

1. AWS Batch User Guide: <https://docs.aws.amazon.com/batch/latest/userguide/what-is-batch.html>
2. AWS DataSync Documentation: <https://aws.amazon.com/datasync/>
3. Amazon S3 Lifecycle Configuration: <https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lifecycle-mgmt.html>
4. Amazon EventBridge User Guide: <https://docs.aws.amazon.com/eventbridge/latest/userguide/eb-what-is.html>
5. AWS Lambda Documentation: <https://docs.aws.amazon.com/lambda/latest/dg/welcome.html>
6. Amazon EC2 Spot Instances: <https://aws.amazon.com/ec2/spot/>
7. AWS VPC Endpoints: <https://docs.aws.amazon.com/vpc/latest/privatelink/vpc-endpoints.html>
8. AWS IAM Best Practices: <https://docs.aws.amazon.com/IAM/latest/UserGuide/best-practices.html>
