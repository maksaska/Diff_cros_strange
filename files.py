import pandas as pd


class Files:
    def __init__(self, name: str):
        # Check the channel
        assert name == "K+L" or name == "K+S0", f"Channel name: {name} is invalid"
        # Data files reading
        if name == "K+L":
            dir = "KL"
        elif name == "K+S0":
            dir = "KS"

        self.__name = name

        Data1 = pd.read_csv(f'./Data/{dir}/P1.csv')
        Data1['St'] = Data1['Su']/(1 + Data1['eps']/5)
        Data1['dSt'] = Data1['dSu']/(1 + Data1['eps']/5)
        Data1['Sl'] = Data1['Su']/(5 + Data1['eps'])
        Data1['dSl'] = Data1['dSu']/(5 + Data1['eps'])
        Data1['E'] = 2.567
        Data1 = Data1.drop(columns=['Su', 'dSu'])

        Data2 = pd.read_csv(f'./Data/{dir}/P2.csv')
        Data2['St'] = Data2['Su']/(1 + Data2['eps']/5)
        Data2['dSt'] = Data2['dSu']/(1 + Data2['eps']/5)
        Data2['Sl'] = Data2['Su']/(5 + Data2['eps'])
        Data2['dSl'] = Data2['dSu']/(5 + Data2['eps'])
        Data2['E'] = 4.056
        Data2 = Data2.drop(columns=['Su', 'dSu'])

        Data3 = pd.read_csv(f'./Data/{dir}/P3.csv')
        Data3['St'] = Data3['Su']/(1 + Data3['eps']/5)
        Data3['dSt'] = Data3['dSu']/(1 + Data3['eps']/5)
        Data3['Sl'] = Data3['Su']/(5 + Data3['eps'])
        Data3['dSl'] = Data3['dSu']/(5 + Data3['eps'])
        Data3['E'] = 5.499
        Data3 = Data3.drop(columns=['Su', 'dSu'])

        self.__Data = pd.concat([Data1, Data2, Data3], ignore_index=True)

        del Data1, Data2, Data3

        self.__Data_Sigma = pd.read_csv(f'./Data/{dir}/Sigma_Photo.csv')
        self.__Data_Diff = pd.read_csv(f'./Data/{dir}/Diff_Photo.csv')
        self.__Data_Q2cos = pd.read_csv(f'./Data/{dir}/K.csv')
        self.__W_syserr = pd.read_csv(f'./Data/{dir}/K_W_syserr.csv')

    @property
    def name(self):
        return self.__name

    @property
    def Data(self):
        return self.__Data

    @property
    def Data_Sigma(self):
        return self.__Data_Sigma

    @property
    def Data_Diff(self):
        return self.__Data_Diff

    @property
    def Data_Q2cos(self):
        return self.__Data_Q2cos

    @property
    def W_syserr(self):
        return self.__W_syserr

    @property
    def cosQ2_trig(self):
        return self.__cosQ2_trig

    def __repr__(self):
        return f'Files("{self.name}")'

    def __str__(self):

        print(self.Data)
        print(self.Data_Sigma)
        print(self.Data_Diff)
        print(self.Data_Q2cos)
        print(self.cosQ2_trig)
        print(self.W_syserr)

        return f'This is all the data for {self.name}'
