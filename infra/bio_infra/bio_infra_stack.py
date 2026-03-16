from aws_cdk import (
    Stack,
    aws_ec2 as ec2,
    aws_batch as batch,
    aws_iam as iam,
)
from constructs import Construct
from typing import cast # 型を明示するために追加

class BioInfraStack(Stack):

    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)

        # 1. VPC の作成
        vpc = ec2.Vpc(
            self, "BioInfraVpc",
            max_azs=2,
        )

        # 2. IAM ロールの作成
        batch_job_role = iam.Role(
            self, "BatchJobRole",
            assumed_by=iam.ServicePrincipal("ecs-tasks.amazonaws.com"),
        )

        bucket_arn = "arn:aws:s3:::ci-genome-tokyo"
        batch_job_role.add_to_policy(
            iam.PolicyStatement(
                actions=["s3:*"],
                resources=[bucket_arn, f"{bucket_arn}/*"],
            )
        )

        # 3. 起動テンプレートの作成 (L1パーツを使用)
        # 変数名を少し変えて、衝突を避けます
        my_lt = ec2.CfnLaunchTemplate(
            self, "BatchLaunchTemplate",
            launch_template_name="bio-batch-storage-lt",
            launch_template_data=ec2.CfnLaunchTemplate.LaunchTemplateDataProperty(
                block_device_mappings=[
                    ec2.CfnLaunchTemplate.BlockDeviceMappingProperty(
                        device_name="/dev/xvda",
                        ebs=ec2.CfnLaunchTemplate.EbsProperty(
                            volume_size=200,      # 目標の 200GB
                            volume_type="gp3"
                        )
                    )
                ]
            )
        )

        compute_env_instance_role = iam.Role(
            self, "ComputeEnvInstanceRole",
            assumed_by=iam.ServicePrincipal("ec2.amazonaws.com"),
            managed_policies=[
                iam.ManagedPolicy.from_aws_managed_policy_name("service-role/AmazonEC2ContainerServiceforEC2Role"),
                iam.ManagedPolicy.from_aws_managed_policy_name("AmazonEC2ContainerRegistryReadOnly"),
            ]
        )

        # 4. 計算環境の作成
        compute_env = batch.ManagedEc2EcsComputeEnvironment(
            self, "BioInfraComputeEnv",
            vpc=vpc,
            spot=True,
            allocation_strategy=batch.AllocationStrategy.SPOT_CAPACITY_OPTIMIZED,
            instance_role=compute_env_instance_role
        )

        # --- 重要：エスケープハッチ (赤い波線を消すための書き方) ---
        # cfn_ce に型を教えることで、エディタの警告を消します
        cfn_ce = cast(batch.CfnComputeEnvironment, compute_env.node.default_child)
        
        # CloudFormationの正式なプロパティ名（大文字開始）で設定を流し込みます
        cfn_ce.add_property_override("ComputeResources.LaunchTemplate", {
            "LaunchTemplateName": "bio-batch-storage-lt",
            "Version": "$Latest"
        })

        # 5. ジョブキューの作成
        batch.JobQueue(
            self, "BioInfraJobQueue",
            job_queue_name="bio-analysis-queue",
            compute_environments=[
                batch.OrderedComputeEnvironment(
                    compute_environment=compute_env,
                    order=1
                )
            ],
        )

