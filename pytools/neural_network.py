import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(X_norm.shape[0], 64),
            # nn.Linear(64, 64),
            # nn.Linear(30, 10),
            nn.Linear(64, 1)
        )

    def forward(self, x):
        x = F.relu(self.model[0](x))
        for layer in self.model[1:]:
            x = F.relu(layer(x))
        return x


