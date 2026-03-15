import aws_cdk as core
import aws_cdk.assertions as assertions

from bio_infra.bio_infra_stack import BioInfraStack

# example tests. To run these tests, uncomment this file along with the example
# resource in bio_infra/bio_infra_stack.py
def test_sqs_queue_created():
    app = core.App()
    stack = BioInfraStack(app, "bio-infra")
    template = assertions.Template.from_stack(stack)

#     template.has_resource_properties("AWS::SQS::Queue", {
#         "VisibilityTimeout": 300
#     })
