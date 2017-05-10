from abc import ABCMeta, abstractmethod
import yaml


class abstractclassmethod(classmethod):

    __isabstractmethod__ = True

    def __init__(self, callable):
        callable.__isabstractmethod__ = True
        super(abstractclassmethod, self).__init__(callable)


class Default(object):
    __metaclass__ = ABCMeta

    @abstractclassmethod
    def from_dict(cls, dct):
        pass

    @abstractmethod
    def to_dict(self):
        pass

    @classmethod
    def import_file(cls, file_name):
        f = open(file_name, "r")
        dict_data = yaml.load(f)
        f.close()
        return cls.from_dict(dict_data)

    def export_file(self, file_name):
        dict_data = self.to_dict()
        f = open(file_name, "w")
        yaml.dump(dict_data, f)
        f.close()
