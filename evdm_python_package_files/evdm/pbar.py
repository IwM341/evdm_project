class progress_bar_jupiter():
    """
    Class for logging progress in jupyter notebook
    """
    def __init__(self,name='Progress'):
        """
        init constructor:
        name -- label
        """
        self.inited = False
        self.name = name
        self.label = None
        self.progress = None
        self.max_count = None
        self.index= None

    def _progress_string(self):
        m_name = self.name
        m_percent = (self.index*100)/self.max_count
        m_time = self._get_time() - self.start_time
        s_prefix = ""
        if (m_time < 1):
            s_prefix = "m"
            m_time = m_time*1000
         
        return f"{m_name} : {m_percent:.2f}%, {self.index}/{self.max_count }, {m_time:.2f}{s_prefix}s"

    def _init_impl(self,max_count):
        from ipywidgets import IntProgress, HTML, VBox
        from IPython.display import display
        self.progress = IntProgress(min=0, max=max_count, value=0)
        self.index = 0
        self.max_count = max_count
        self.inited = True
        try:
            import time
            self._get_time = time.time
        except:
            self._get_time = lambda: 0
        self.start_time = self._get_time()

        self.label = HTML(self._progress_string())
        self.box = VBox(children=[self.label, self.progress])
        display(self.box)
    
    def update(self,current_progress,full_progress):
        """
        function, which updates current state
        """
        if(self.inited):
            self.index = current_progress
            self.label.value = self._progress_string()
            self.progress.value = self.index
        else:
            self._init_impl(full_progress)
class progress_bar_cmd():
    """
    Class for logging progress by simple print
    """
    def __init__(self,name='Progress'):
        """
        init constructor:
        name -- label
        """
        self.inited = False
        self.name = name
        self.max_count = None
        self.index= None
       
    
    def _print_progress(self):
        m_name = self.name
        m_percent = (self.index*100)/self.max_count
        print(f'\r{m_name}: {m_percent}%, {self.index}/{self.max_count}')

    def _init_impl(self,max_count):
        self.max_count = max_count
        self.index = 0
        self.inited = True
        self._print_progress()
        
    
    def update(self,current_progress,full_progress):
        """
        function, which updates current state
        """
        if(self.inited):
            self.index = current_progress
            self._print_progress()
        else:
            self._init_impl(full_progress)