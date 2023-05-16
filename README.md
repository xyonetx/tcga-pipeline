# TCGA processing pipeline

## To run:

**Prerequisites:**
- AWS account
- Nextflow is installed (https://www.nextflow.io/docs/latest/getstarted.html)

Many components cannot be run on local hardware. We expect that you will be executing this pipeline tool using an AWS Batch runner. Therefore, you must first set up an AWS Batch queue using the repository here: https://github.com/xyonetx/batch (or via some other method).

Once the Batch environment is ready, you must configure Nextflow. Copy `nextflow.batch.config.tmpl` and fill in the variables. Note:
- `process.queue` is the name of the Batch queue you created
- `aws.batch.cliPath` should be changed if your AWS cli tool is *not* located at `/opt/aws-cli/bin/aws`. If you used https://github.com/xyonetx/batch to configure your Batch queue, you should have made an EC2 AMI which included this tool, ideally at this location. Modify the path as necessary, taking care to note issues related to shadowing of Docker paths if the AWS cli is installed under `/usr/bin` or `/usr/local/bin`.
- `workDir` should be an S3 bucket (or folder in bucket) where AWS Batch can write. As part of setting up Batch, you created a bucket- just use that one. Include the `s3://` prefix.

Next we set the pipeline parameters. These might change over time, so you can consult the meaning of most parameters by scanning the `main.tf` script. Copy `params.json.tmpl` and fill in variables as required.

Now you can run with:
```
nextflow run main.nf -c <NEXTFLOW CONFIG> -params-file <PARAMS JSON>
```