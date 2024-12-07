"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

QtWidgets.py
Extensions and additions to PyQt5.QtWidgets for the RxGaming tool.

This is a collection of generally unrelated widgets and widget helpers. They all extend PyQt
in some way and act like Qt elements, so they are here together.
"""

import math
import random
import string
from enum import Enum, auto
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QPushButton, QFileDialog, QSlider, QMainWindow
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QPoint
from PyQt5.QtGui import QPainter, QBrush, QColor, QPen, QFont, QFontMetrics


class QFileSelectionLineEdit(QWidget):
    """Line edit together with browse button to select files or directories."""

    class FileType(Enum):
        """Enumeration for indicating whether to open a file or a directory."""
        
        File = auto()
        """Open a file."""
        Directory = auto()
        """Open a directory."""
    
    def __init__(self, caption="Browse...", file_type=FileType.File, filter="Any file (*.*)", new_file=False):
        """Create a new ``QFileSelectionLineEdit``
        
        The file selection dialog window title is set by ``caption``. Pass in ``file_type`` to
        select whether to browse for files or directories. When browsing for files, ``filter``
        can be set to a standard filter text string (or list of strings) to limit possible choices.
        """
        
        super().__init__()
        self.caption = caption
        self.file_type = file_type
        self.filter = filter
        self.new_file = new_file
        
        self.layout = QHBoxLayout()
        self.line_edit = QLineEdit()
        self.button = QPushButton("Browse")
        self.layout.addWidget(self.line_edit)
        self.layout.addWidget(self.button)
        self.setLayout(self.layout)
        
        self.button.clicked.connect(self._on_browse)
    
    def _on_browse(self, button):
        # TODO: start from last place ended up
        if self.file_type is self.FileType.File:
            if self.new_file:
                file_path = QFileDialog.getSaveFileName(self, self.caption, "", self.filter)[0]
            else:
                file_path = QFileDialog.getOpenFileName(self, self.caption, "", self.filter)[0]
        else:
            file_path = QFileDialog.getExistingDirectory(self, self.caption, "")
        if not file_path == "":
            self.line_edit.setText(file_path)
    
    def text(self):
        """Get text from the line edit."""
        
        return self.line_edit.text()
    
    def setText(self, text):
        """Set text in the line edit."""
        self.line_edit.setText(text)


class QWindow(QWidget):
    """Extends QWidget for the purpose of hooking in to the window closure callback.
    
    To be used exactly like a QWidget that is a window (has no parent).
    """
    
    def __init__(self):
        super().__init__()
        self.onClosed = None
    
    def setOnClosed(self, callback):
        """Set callback for window close event.
        
        The ``callback`` should accept one argument, the ``QWindow`` that is closing.
        """
        
        self.onClosed = callback
    
    def closeEvent(self, event):
        """"""
        super().closeEvent(event)
        if self.onClosed is not None:
            self.onClosed(self)


class QMainWindowRx(QMainWindow):
    def __init__(self):
        super().__init__()
        self.onClosed = None

    def setOnClosed(self, callback):

        self.onClosed = callback

    def closeEvent(self, event):
        super().closeEvent(event)
        if self.onClosed is not None:
            self.onClosed(self)


class QWaitingIndicator(QWidget):
    """Add waiting spinny thing overlay to a widget.
    
    This widget is hidden by default. When ``show()`` is called, it appears
    on top of its parent widget with a little waiting circle animation thing.
    Optionally, text can be displayed above the animation.
    
    This widget is best used as an overlay on another widget, for example::
    
        window = QWindow()
        ...  # set up window elements
        start_some_processing_on_another_thread()
        window.overlay = QWaitingIndicator(window)
        window.overlay.setDisplayText("Please wait...")
        window.overlay.show()
        
        ...  # elsewhere, once processing is done
        window.overlay.hide()
    """
    
    def __init__(self, parent):
        super().__init__(parent)
        parent.resizeEvent = self.parentResizeEvent
        self.displayText = ""
        self.counter = 0
        self.hide()
    
    def setDisplayText(self, displayText):
        """Set text to be displayed above animation"""
        
        self.displayText = displayText
    
    def parentResizeEvent(self, event):
        self.resize(event.size())
        event.accept()
    
    def _blendColor(self, value1, value2, ratio):
        """Single-band color blending utility function.
        
        When ``ratio = 0``, return is ``value1``.
        When ``ratio = 1``, return is ``value2``.
        Values between 0 and 1 blend correspondingly.
        """
        
        return value1 + int(ratio * float(value2 - value1))
    
    def _blendColors(self, color1, color2, ratio=0.5):
        """RGB color blending utility function."""
        
        return QColor(self._blendColor(color1.red(), color2.red(), ratio),
                      self._blendColor(color1.green(), color2.green(), ratio),
                      self._blendColor(color1.blue(), color2.blue(), ratio),
                      self._blendColor(color1.alpha(), color2.alpha(), ratio))
    
    def paintEvent(self, event):
        """"""
        painter = QPainter()
        painter.begin(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.fillRect(event.rect(), QBrush(QColor(255, 255, 255, 127)))
        painter.setPen(QPen(Qt.NoPen))
        
        # TODO make these properties
        n_dots = 12
        circle_dia = 25
        dot_dia = 10
        tail_length = 4
        lit_col = QColor(255, 100, 100)
        dull_col = QColor(100, 100, 100)
        
        colored_dots = [(self.counter - i) % n_dots for i in range(tail_length + 1)]
        ratio_step = 1 / (tail_length + 2)
        ratios = {colored_dots[i]: ratio_step * i for i in range(tail_length + 1)}
        for i in range(n_dots):
            if i in colored_dots:
                painter.setBrush(QBrush(self._blendColors(lit_col, dull_col, ratios[i])))
            else:
                painter.setBrush(QBrush(dull_col))
            painter.drawEllipse(
                self.width() / 2 + circle_dia * math.cos(2 * math.pi * i / float(n_dots)),
                self.height() / 2 + circle_dia * math.sin(2 * math.pi * i / float(n_dots)),
                dot_dia, dot_dia)
        
        painter.setPen(QColor(30, 30, 30))  # TODO make font settings a property
        font = painter.font()
        font.setFamily("Helvetica")
        font.setPointSize(12)
        font.setWeight(QFont.Bold)
        painter.setFont(font)
        rect = painter.fontMetrics().boundingRect(self.displayText)
        rect.moveTo(int((self.width() - rect.width()) / 2),
                    int((self.height() - rect.height()) / 2 - circle_dia - 2 * dot_dia))
        painter.drawText(rect, Qt.AlignCenter, self.displayText)
        painter.end()
    
    def timerEvent(self, event):
        """"""
        self.counter += 1
        self.update()
    
    def showEvent(self, event):
        """"""
        self.counter = 0
        self.timer = self.startTimer(150)  # TODO make speed a property
    
    def hideEvent(self, event):
        """"""
        self.killTimer(self.timer)


class SliderWithValue(QSlider):
    def __init__(self, parent=None):
        super(SliderWithValue, self).__init__(parent)

        self.stylesheet = """
        QSlider::groove:vertical {
                background-color: #222;
                width: 30px;
        }
        QSlider::handle:vertical {
            border: 1px #438f99;
            border-style: outset;
            margin: -2px 0;
            width: 30px;
            height: 3px;
            background-color: #438f99;
        }
        QSlider::sub-page:vertical {
            background: #4B4B4B;
        }
        QSlider::groove:horizontal {
                background-color: #222;
                height: 30px;
        }
        QSlider::handle:horizontal {
            border: 1px #438f99;
            border-style: outset;
            margin: -2px 0;
            width: 3px;
            height: 30px;
            background-color: #438f99;
        }
        QSlider::sub-page:horizontal {
            background: #4B4B4B;
        }
        """

        # self.setStyleSheet(self.stylesheet)

    def paintEvent(self, event):
        QSlider.paintEvent(self, event)

        curr_value = str(self.value())
        round_value = round(float(curr_value), 2)

        painter = QPainter(self)
        painter.setPen(QPen(Qt.black))

        font_metrics = QFontMetrics(self.font())
        font_width = font_metrics.boundingRect(str(round_value)).width()
        font_height = font_metrics.boundingRect(str(round_value)).height()

        rect = self.geometry()
        if self.orientation() == Qt.Horizontal:
            horizontal_x_pos = rect.width() - font_width - 5
            horizontal_y_pos = rect.height() * 0.75

            painter.drawText(QPoint(horizontal_x_pos, horizontal_y_pos), str(round_value))

        elif self.orientation() == Qt.Vertical:
            vertical_x_pos = rect.width() - font_width - 5
            vertical_y_pos = rect.height() * 0.75

            painter.drawText(QPoint(rect.width() / 2.0 - font_width / 2.0, rect.height() - 5), str(round_value))
        else:
            pass

        painter.drawRect(rect)


def pyqtWorkerthread(return_type=type(None), slot=None):
    """Decorator, moves processing off the UI thread.
    
    The decorated method will be spun up on its own thread when called.
    It will immediately return the ``QThread`` object associated with the new thread.
    To get the returned value of the decorated function itself, it is necessary
    to designate a slot that will be called upon completion.
    PyQt requires that the return type of a signal is well-known,
    so ``return_type`` must be set to the type of object returned by the decorated method
    for the return slot to work.
    The only supported paradigm is for ``slot`` to be another method within
    the same class as the decorated method, and declared before the decorated method.
    
    .. note:: Use this decorator only on instance methods. The return slot must be an instance
       method as well.
    
    Here is an example of a worker thread with no return::
    
        def someMethod(self, piece_of_work):
            self.do_work(piece_of_work)
        
        @pyqtWorkerthread()
        def do_work(self, material)
            ...  # do something with material
    
    And an example with a return and UI update::
    
        def someMethod(self, piece_of_work):
            self.do_work(piece_of_work)
        
        def work_done(self, product_value):
            self.label.setText(str(product_value))
        
        @pyqtWorkerthread(float, work_done)
        def do_work(self, material):
            ...  # do something with material
            return float(material.result)
    """
    
    def decorator_workerthread(func):
        def dowork(self, *args, **kwargs):
            thread_name = "_thread" + "".join(random.choices(string.ascii_letters + string.digits, k=16))

            class WorkerThread(QThread):
                name = thread_name
                finished_signal = pyqtSignal(return_type)

                def run(self):
                    self.finished_signal.emit(func(self, *args, **kwargs))
                    delattr(self.parent, thread_name)
            worker_instance = WorkerThread()
            worker_instance.parent = self
            setattr(self, thread_name, worker_instance)
            if slot is not None:
                worker_instance.finished_signal.connect(getattr(self, slot.__name__))
            worker_instance.start()
            return worker_instance
        return dowork
    return decorator_workerthread
